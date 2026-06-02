# pacotes
BiocManager::install("edgeR")
install.packages("openxlsx")
install.packages("ggrepel")
install.packages("ggridges")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

# Carrega a biblioteca para ler arquivos do Excel
library(readxl)
library(edgeR)
library(limma)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db) # Banco de dados do Genoma Humano
library(enrichplot)
library(ggplot2)
library(tibble)

# (Opcional) Define a pasta de trabalho onde estão seus arquivos
setwd("C:/Users/conta/OneDrive/doutorado/analises_exp_rnaseq")

# 1. Importa a matriz de ABA (Aba 1)
expressao <- read_excel("RNAseq_results_HMOX1KO_mock_and_addback_genes.xlsx", sheet = 1)

# 2. Importa o metadata (Aba 2)
metadata <- read_excel("raw_clarissa.xlsx", sheet = 2)



# 2. Define as colunas de anotação que você quer manter
# (Ajuste os nomes conforme aparecem exatamente no seu Excel)
cols_info <- c("Chr", "gene", "gene_biotype", "GeneID", "description")

# 3. Define as colunas de expressão (as suas amostras)
# Você pode listar manualmente ou usar um padrão (ex: que começam com A375)
todas_amostras <- c(
  "A375_HMOX1_addback_Ctrl_1",  "A375_HMOX1_addback_Ctrl_2",  "A375_HMOX1_addback_Ctrl_3",
  "A375_HMOX1_addback_Hemin_1", "A375_HMOX1_addback_Hemin_2", "A375_HMOX1_addback_Hemin_3",
  "A375_HMOX1_mock_Ctrl_1",     "A375_HMOX1_mock_Ctrl_2",     "A375_HMOX1_mock_Ctrl_3",
  "A375_HMOX1_mock_Hemin_1",    "A375_HMOX1_mock_Hemin_2",    "A375_HMOX1_mock_Hemin_3"
)

# Amostras descartadas devido a problemas de qualidade biológica/técnica
outliers_identificados <- c("A375_HMOX1_addback_Ctrl_1", "A375_HMOX1_addback_Hemin_3")

# cols_expr conterá apenas as 10 amostras limpas e válidas para a estatística
cols_expr <- setdiff(todas_amostras, outliers_identificados)

# 4. Criamos a matriz numérica apenas com as amostras selecionadas
matriz_num <- as.matrix(expressao[, cols_expr])
rownames(matriz_num) <- make.unique(as.character(expressao$gene))

# 5. Guardamos as informações dos genes em um dataframe separado para o final
genes_info <- expressao[, cols_info]

# Garantir que o metadata só tenha as amostras que você selecionou em 'cols_expr'
metadata <- metadata %>%
  filter(SampleID %in% cols_expr) %>%
  arrange(match(SampleID, colnames(matriz_num))) # Garante a mesma ordem

# Verificação de segurança
if(all(colnames(matriz_num) == metadata$SampleID)) {
  print("✅ Sincronização OK!")
} else {
  stop("❌ Erro: As colunas da matriz não batem com as linhas do metadata.")
}


## LIMPEZA DA TABELA: EXCLUSÃO DE VALORES NULOS ################################

# 1. Definição do corte (Threshold)
# Usamos 1 pois os dados são TPM/CPM (escala linear)
threshold <- 5

# 2. Regra do menor grupo
# Se você tem 3 replicatas por grupo, o gene deve estar expresso em pelo menos 3 amostras
min_samples <- 2 

# 3. Filtragem
keep <- rowSums(matriz_num >= threshold) >= min_samples # se expresso em no mínimo 2 amostras
matriz_filtrada <- matriz_num[keep, ]

# 4. LOG-TRANSFORMAÇÃO (Obrigatório para limma)
# Adicionamos +1 para evitar o erro matemático de Log de zero (log2(0) é impossível)
matriz_log <- log2(matriz_filtrada + 1)

# Verificação
print(paste("Original:", nrow(matriz_num), "genes"))
print(paste("Após filtragem:", nrow(matriz_filtrada), "genes"))

### matriz com outliers (para pca) ################################################
matriz_com_outliers <- as.matrix(expressao[, todas_amostras])
rownames(matriz_com_outliers) <- make.unique(as.character(expressao$gene))

# 2. Aplicar a mesma filtragem e log-transformação básica para consistência
keep_all <- rowSums(matriz_com_outliers >= threshold) >= 3
matriz_filtrada_all <- matriz_com_outliers[keep_all, ]
matriz_log_all <- log2(matriz_filtrada_all + 1)

# 3. Remover genes com variância zero usando funções nativas do R
variancia_all <- apply(matriz_log_all, 1, var)
matriz_log_all_filtrada <- matriz_log_all[variancia_all > 0, ]

###GRÁFICO DE PCA ##############################################################

cores_pca <- c(
  "Addback_NT"  = "#8A2BE2", # Roxo BlueViolet
  "Addback_Hemin" = "#4B0082", # Roxo Indigo (Mais escuro)
  "Mock_NT"     = "#4169E1", # Azul Real
  "Mock_Hemin"    = "#000080"  # Azul Marinho Escuro
)

pca_all <- prcomp(t(matriz_log_all_filtrada), scale. = TRUE)

# 5. Montar a tabela de dados para o ggplot e criar a coluna de Grupo para a cor
df_pca_all <- as.data.frame(pca_all$x) %>%
  rownames_to_column("SampleID") %>%
  left_join(read_excel("raw_clarissa.xlsx", sheet = 2), by = "SampleID") %>%
  mutate(
    Is_Outlier = ifelse(SampleID %in% outliers_identificados, "Outlier", "Válida"),
    Grupo_Cor = paste(Genotype, Treatment, sep = "_")
  )

# Calcular a variância explicada
var_all <- round(100 * pca_all$sdev^2 / sum(pca_all$sdev^2), 1)

# 6. Gráfico 1: Com Outliers
plot_pca_com_outliers <- ggplot(df_pca_all, aes(x = PC1, y = PC2, color = Grupo_Cor, shape = Grupo_Cor)) +
  geom_point(size = 4, alpha = 0.8) + # Aqui todos os pontos ganham tamanho 4 de forma justa!
  geom_text_repel(aes(label = SampleID), size = 3, show.legend = FALSE) +
  scale_color_manual(values = cores_pca, name = "Grupo / Tratamento") +
  scale_shape_manual(values = c("Addback_NT" = 17, "Addback_Hemin" = 16, "Mock_NT" = 17, "Mock_Hemin" = 16), name = "Grupo / Tratamento") +
  labs(
    title = "Principal Component Analysis (PCA) of Full Transcriptomic Profiles",
    subtitle = "Identification of sample outliers based on global variance",
    x = paste0("PC1 (", var_all[1], "%)"),
    y = paste0("PC2 (", var_all[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "right"
  )

### PCA 2: REMOVENDO OS OUTLIERS (DADOS LIMPOS) --------------------------------

# 1. Remover genes com variância zero na sua matriz_log limpa
variancia_clean <- apply(matriz_log, 1, var)
matriz_log_clean_filtrada <- matriz_log[variancia_clean > 0, ]

# 2. Rodar o PCA na matriz limpa
pca_clean <- prcomp(t(matriz_log_clean_filtrada), scale. = TRUE)

# 3. Montar a tabela de dados e criar a coluna de Grupo para a cor/forma
df_pca_clean <- as.data.frame(pca_clean$x) %>%
  rownames_to_column("SampleID") %>%
  left_join(metadata, by = "SampleID") %>%
  mutate(Grupo_Cor = paste(Genotype, Treatment, sep = "_"))

# Calcular a variância explicada
var_clean <- round(100 * pca_clean$sdev^2 / sum(pca_clean$sdev^2), 1)

# 4. Gráfico 2: Sem Outliers (Versão Unificada com Legenda Limpa)
plot_pca_sem_outliers <- ggplot(df_pca_clean, aes(x = PC1, y = PC2, color = Grupo_Cor, shape = Grupo_Cor)) +
  geom_point(size = 4, alpha = 0.8) + # Tamanho 4 fixo para combinar com o Gráfico A
  geom_text_repel(aes(label = SampleID), size = 3, show.legend = FALSE) +
  scale_color_manual(values = cores_pca, name = "Grupo / Tratamento") +
  # 17 = Triângulo (NT), 16 = Círculo (Hemin)
  scale_shape_manual(values = c("Addback_NT" = 17, "Addback_Hemin" = 16, "Mock_NT" = 17, "Mock_Hemin" = 16), name = "Grupo / Tratamento") +
  labs(
    title = "Principal Component Analysis (PCA) of Filtered RNA-Seq Dataset",
    subtitle = "Distinct clustering patterns by genotype and treatment after outlier removal",
    x = paste0("PC1 (", var_clean[1], "%)"),
    y = paste0("PC2 (", var_clean[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "right"
  )


### EXIBIR OS GRÁFICOS --------------------------------------------------------
plot(plot_pca_com_outliers)
plot(plot_pca_sem_outliers)

## ANÁLISE MATRIZ DE DESIGN ####################################################

# 1. Criar um fator que combina Genótipo e Tratamento
# Isso vai gerar etiquetas como: Mock_NT, Mock_Hemin, Addback_NT, Addback_Hemin
grupo <- factor(paste(metadata$Genotype, metadata$Treatment, sep = "_"))

# 2. Criar a Design Matrix
# O "~ 0" diz ao R para não assumir um 'padrão' (intercepto) e tratar todos os grupos de forma justa.
design <- model.matrix(~ 0 + grupo)

# Ajustar os nomes das colunas para facilitar a leitura
colnames(design) <- levels(grupo)

# Ver a matriz (opcional)
print(design)

## ANÁLISE: LIMMA FIT ##########################################################

# Ajusta o modelo linear aos dados logaritmizados
fit <- lmFit(matriz_log, design)


## ANÁLISES: DESENHO EXPERIMENTAL ##############################################

# Definimos os contrastes (Comparações)
cont.matrix <- makeContrasts(
  # Pergunta A: Qual o papel da HMOX1 sem nenhum tratamento?
  HMOX1_Basal = Addback_NT - Mock_NT,
  
  # Pergunta B: O que o Hemin faz no Mock (independente de HMOX1)?
  Hemin_Effect_Mock = Mock_Hemin - Mock_NT,
  
  # Pergunta C: O que o Hemin faz no Add-back (onde a HMOX1 está presente)?
  Hemin_Effect_Addback = Addback_Hemin - Addback_NT,
  
  # Pergunta D (A INTERAÇÃO): Onde a resposta ao Hemin é DIFERENTE por causa da HMOX1?
  Interaction = (Addback_Hemin - Addback_NT) - (Mock_Hemin - Mock_NT),
  
  treated_addback_x_mock = Mock_Hemin - Addback_Hemin,
  
  levels = design
)

# Aplicar os contrastes ao modelo
fit2 <- contrasts.fit(fit, cont.matrix)

# ATUALIZAÇÃO: Estatística robusta com correção de tendência
# trend = TRUE: modela a relação entre a média e a variância, essencial para RNA-seq
# robust = TRUE: reduz o peso de outliers sem precisar excluir amostras do seu estudo
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

## 8. EXTRAÇÃO E FILTRAGEM DOS RESULTADOS SIGNIFICATIVOS #####################################

# 1. Preparar a tabela de anotação para o "Join"
# Criamos uma chave única (gene_ID_unico) para bater exatamente com os rownames da matriz
genes_info_clean <- genes_info %>%
  mutate(gene_ID_unico = make.unique(as.character(gene)))

# 2. Função para extrair, anotar e filtrar resultados
extrair_e_anotar <- function(fit_obj, nome_coef, info_ref) {
  topTable(fit_obj, coef = nome_coef, number = Inf, sort.by = "P") %>%
    rownames_to_column("gene_ID_unico") %>%
    left_join(info_ref, by = "gene_ID_unico") %>%
    # Filtro de Qualidade: Proteínas, Significância e Fold Change
    filter(gene_biotype == "protein_coding" & adj.P.Val < 0.05 & abs(logFC) > 1) %>%
    arrange(adj.P.Val)
}

# 3. Gerar as listas de genes significativos
sig_basal     <- extrair_e_anotar(fit2, "HMOX1_Basal", genes_info_clean)
sig_mock      <- extrair_e_anotar(fit2, "Hemin_Effect_Mock", genes_info_clean)
sig_addback   <- extrair_e_anotar(fit2, "Hemin_Effect_Addback", genes_info_clean)
sig_interacao <- extrair_e_anotar(fit2, "Interaction", genes_info_clean)
sig_tratadas <- extrair_e_anotar(fit2, "treated_addback_x_mock", genes_info_clean)



## 9. EXPORTAÇÃO PARA EXCEL ###################################################

# Criar a lista de abas
lista_final <- list(
  "Interacao_HMOX1_Hemin" = sig_interacao,
  "Efeito_Basal_HMOX1"    = sig_basal,
  "Hemin_no_Mock"         = sig_mock,
  "Hemin_no_Addback"      = sig_addback
)

# Salvar arquivo (Altere o nome se desejar)
write.xlsx(lista_final, file = "Resultados_Sincronizados_RNAseq_Bela.xlsx")

print("✅ Arquivo Excel com resultados significativos gerado com sucesso!")

## 9. LOGICA PARA GRÁFICOS (VOLCANO E BARRAS) ##################################

# 1. Função para extrair dados COMPLETOS (Necessário para o fundo do Volcano)
extrair_full <- function(fit_obj, nome_coef, info_ref) {
  topTable(fit_obj, coef = nome_coef, number = Inf, sort.by = "P") %>%
    rownames_to_column("gene_ID_unico") %>%
    left_join(info_ref, by = "gene_ID_unico") %>%
    filter(gene_biotype == "protein_coding") %>%
    mutate(Status = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE ~ "Non-significant"
    ))
}

# 2. Gerar tabelas completas para os gráficos
full_basal     <- extrair_full(fit2, "HMOX1_Basal", genes_info_clean) 
full_hemin_mock <- extrair_full(fit2, "Hemin_Effect_Mock", genes_info_clean)
full_hemin_addback <- extrair_full(fit2, "Hemin_Effect_Addback", genes_info_clean)
full_hemin_both <- extrair_full(fit2, "treated_addback_x_mock", genes_info_clean)
full_interacao <- extrair_full(fit2, "Interaction", genes_info_clean) %>%
  mutate(Status = case_when(
    adj.P.Val < 0.05 & logFC > 1 ~ "Positive Interaction",
    adj.P.Val < 0.05 & logFC < -1 ~ "Negative Interaction",
    TRUE ~ "Non-significant"
  )) 

## FUNÇÃO DIVA PARA VOLCANO PLOTS #############################################

gerar_volcano_layout_final <- function(dados_completos, dados_sig, titulo = "", subtitulo = "DEG | padj < 0.05 | |log2FC| > 1", n_cada_lado = 10) {
  
  # 1. Seleciona o Top N de cada lado (Up e Down) com base no p-valor
  top_up <- dados_sig %>% 
    filter(logFC > 0) %>% 
    arrange(adj.P.Val) %>% 
    head( n_cada_lado)
  
  top_down <- dados_sig %>% 
    filter(logFC < 0) %>% 
    arrange(adj.P.Val) %>% 
    head( n_cada_lado)
  
  # Une os dois grupos para as etiquetas
  top_labels <- bind_rows(top_up, top_down)
  
  # 2. Contagem total de variáveis
  total_vars <- nrow(dados_completos)
  
  # 3. Construção do Gráfico
  ggplot(dados_completos, aes(x = logFC, y = -log10(adj.P.Val), color = Status)) +
    theme_bw() +
    geom_point(alpha = 1, size = 1.5) +
    
    # Cores conforme o layout solicitado[cite: 1, 2]
    scale_color_manual(values = c(
      "Downregulated" = "#0000FF",        # Azul
      "Negative Interaction" = "#0000FF", 
      "Upregulated" = "#FF0000",          # Vermelho
      "Positive Interaction" = "#FF0000",
      "Non-significant" = "#BEBEBE"       # Cinza
    )) +
    
    # Linhas de Threshold[cite: 1, 2]
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.4) +
    
    # Labels em caixinhas com setas equilibradas[cite: 1, 2]
    geom_label_repel(
      data = top_labels,
      aes(label = gene),
      size = 3,
      fontface = "plain",
      color = "black",
      fill = "white",
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey30",
      segment.size = 0.4,
      arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
      max.overlaps = Inf
    ) +
    
    labs(
      title = titulo,
      subtitle = "DEG | padj < 0.05 | |log2FC| > 1",
      x = "log2 Fold Change",
      y = "-log10 (adj. p-value)",
      caption = paste0("total = ", total_vars, " variables")
    ) +
    
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 11, hjust = 0),
      plot.caption = element_text(size = 10, hjust = 1),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11, color = "black"),
      panel.grid.minor = element_line(color = "#f0f0f0"),
      panel.grid.major = element_line(color = "#e0e0e0")
    )
}
## CHAMADA DO GRÁFICO (Exemplo para Interação) #################################

plot_final_interacao <- gerar_volcano_layout_final(
  dados_completos = full_interacao, 
  dados_sig = sig_interacao,
  titulo = "HMOX1 x Hemin Interaction Analysis",
  n_cada_lado = 5 # Ajuste quantos nomes quer que apareçam
)

plot_final_basal <- gerar_volcano_layout_final(
  dados_completos = full_basal, 
  dados_sig = sig_basal,
  titulo = "HMOX1 Basal Effect (Addback x Mock)",
  n_cada_lado = 8 # Ajuste quantos nomes quer que apareçam
)

plot_final_hemin_mock <- gerar_volcano_layout_final(
  dados_completos = full_hemin_mock, 
  dados_sig = sig_mock,
  titulo = "Hemin Response in Mock Cells (Hemin x NT)",
  n_cada_lado = 5 # Ajuste quantos nomes quer que apareçam
)

plot_final_hemin_addback <- gerar_volcano_layout_final(
  dados_completos = full_hemin_addback, 
  dados_sig = sig_addback,
  titulo = "Hemin Response in Addback Cells (Hemin x NT)",
  n_cada_lado = 5 # Ajuste quantos nomes quer que apareçam
)

plot_final_hemin_both <- gerar_volcano_layout_final(
  dados_completos = full_hemin_both, 
  dados_sig = sig_tratadas,
  titulo = "Hemin Response: A375 mock x A375 Add back",
  n_cada_lado = 8 # Ajuste quantos nomes quer que apareçam
)


plot(plot_final_basal)
plot(plot_final_interacao)
plot(plot_final_hemin_mock)
plot(plot_final_hemin_addback)
plot(plot_final_hemin_both)


###   GRÁFICOS DE REDES DE VIAS - ANÁLISE SECUNDÁRIA ###############################
# 1. Preparar a lista ranqueada (O motor do GSEA)
# Pegamos a tabela completa, removemos qualquer NA por segurança e ordenamos do maior logFC pro menor
lista_gsea_basal <- full_basal %>%
  filter(!is.na(logFC)) %>%
  arrange(desc(logFC))

# Criar o vetor numérico nomeado (exigência do pacote)
vetor_gsea_basal <- lista_gsea_basal$logFC
names(vetor_gsea_basal) <- lista_gsea_basal$gene

# 2. Rodar o GSEA focado em Processos Biológicos (GO)
set.seed(123) # Garante que a análise seja reprodutível
gsea_basal <- gseGO(
  geneList     = vetor_gsea_basal,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",      # Focando em Processos Biológicos
  minGSSize    = 10,        # Ignora vias com menos de 10 genes (muito pequenas)
  maxGSSize    = 500,       # Ignora vias com mais de 500 genes (muito genéricas)
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose      = FALSE
)

# Verificação rápida
print(paste("Vias enriquecidas encontradas:", nrow(as.data.frame(gsea_basal))))

grafico_dotplot <- dotplot(gsea_basal, showCategory = 15, font.size = 11) +
  labs(
    title = "Enriched Biological Processes",
    subtitle = "Basal Effect of HMOX1 (Addback vs Mock)"
  ) +
  # Ajuste da paleta de cores para um visual mais suave e limpo
  scale_color_gradient(
    low = "#7DA1C4",  # Azul pastel (mais significativo)
    high = "#D67A7A", # Vermelho pastel (menos significativo)
    name = "p.adjust"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(size = 11, color = "black"),
    text = element_text(family = "sans") # Mantém a tipografia limpa
  )

# Exibir o gráfico
plot(grafico_dotplot)

# Gerar o Dotplot Dividido (Ativados vs Reprimidos) ####################################################
grafico_dotplot_dividido <- dotplot(gsea_basal, showCategory = 10, split = ".sign", font.size = 11) +
  facet_grid(~.sign) + # Isso cria as duas colunas baseadas no sinal do NES
  labs(
    title = "Enriched Biological Processes (Activated vs Suppressed)",
    subtitle = "Basal Effect of HMOX1 (Addback vs Mock)"
  ) +
  # Mantendo a paleta de cores suave
  scale_color_gradient(
    low = "#7DA1C4",  # Azul pastel (mais significativo)
    high = "#D67A7A", # Vermelho pastel (menos significativo)
    name = "p.adjust"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"), # Deixa o título da coluna bem visível
    strip.background = element_rect(fill = "#f0f0f0", color = "black"), # Fundo cinza lofi pro título da coluna
    text = element_text(family = "sans")
  )

# Exibir o gráfico
plot(grafico_dotplot_dividido)


# Instale o pacote caso ainda não tenha rodando: install.packages("pheatmap")
install.packages("pheatmap")
library(pheatmap)

# 1. Selecionar apenas as amostras basais (NT = Ctrl)
# OBS: O Addback_Ctrl_1 já está de fora respeitando a sua limpeza de outliers!
amostras_basais <- c(
  "A375_HMOX1_mock_Ctrl_1", "A375_HMOX1_mock_Ctrl_2", "A375_HMOX1_mock_Ctrl_3",
  "A375_HMOX1_addback_Ctrl_2", "A375_HMOX1_addback_Ctrl_3"
)

# 2. Extrair os Top 30 genes mais significativos da sua tabela sig_basal
top_30_genes <- sig_basal %>%
  arrange(adj.P.Val) %>%
  head(30) %>%
  pull(gene_ID_unico) # Puxa exatamente o ID que está na matriz_log

# 3. Filtrar a matriz de expressão logaritmizada
matriz_heatmap <- matriz_log[top_30_genes, amostras_basais]

# 4. Criar a "tarja" de anotação para o topo do gráfico ficar didático
# 4. Criar a "tarja" de anotação para o topo do gráfico ficar didático
anotacao_colunas <- metadata %>%
  filter(SampleID %in% amostras_basais) %>%
  dplyr::select(SampleID, Genotype) %>%    # O segredo está aqui: dplyr::select
  column_to_rownames("SampleID")

# Cores da tarja de anotação (mantendo os tons do PCA)
cores_anotacao <- list(
  Genotype = c("Mock" = "#4169E1", "Addback" = "#8A2BE2") 
)

# 5. Gerar o Heatmap
grafico_heatmap <- pheatmap(
  matriz_heatmap,
  scale = "row",              # Transforma em Z-score (essencial para comparar genes de alta/baixa exp)
  annotation_col = anotacao_colunas,
  annotation_colors = cores_anotacao,
  color = colorRampPalette(c("#7DA1C4", "#F5F5F5", "#D67A7A"))(50), # Azul -> Branco -> Rosa suave
  cluster_cols = FALSE,       # Não agrupa colunas (Mantém Mock fixo de um lado e Addback do outro)
  show_rownames = TRUE,       # Mostra o nome dos genes
  show_colnames = FALSE,      # Esconde os nomes longos das amostras para ficar mais limpo
  fontsize_row = 10,          # Tamanho da fonte dos genes
  fontsize = 11,
  main = "Top 30 DEGs: Basal Effect (Addback vs Mock)",
  fontfamily = "sans"
)

# =========================================================================
# 10. GSEA KEGG: SEPARANDO UP E DOWN NO KO (MOCK + HEMIN vs ADDBACK + HEMIN)
# =========================================================================

# 1. Preparar a lista ranqueada usando o contraste direto (Mock_Hemin vs Addback_Hemin)
# Valores positivos de logFC = UP no KO (Mock)
# Valores negativos de logFC = DOWN no KO (Mock)
lista_kegg_tratadas <- full_hemin_both %>%
  filter(!is.na(logFC)) %>%
  arrange(desc(logFC))

# Traduzir para Entrez ID (Exigência do banco KEGG)
traducao_kegg_tratadas <- bitr(
  lista_kegg_tratadas$gene, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)

# Unir tradução com logFC e remover duplicatas
df_kegg_tratadas_ready <- lista_kegg_tratadas %>%
  inner_join(traducao_kegg_tratadas, by = c("gene" = "SYMBOL")) %>%
  distinct(ENTREZID, .keep_all = TRUE) %>% 
  arrange(desc(logFC))

# Criar o vetor ranqueado final
vetor_gsea_tratadas <- df_kegg_tratadas_ready$logFC
names(vetor_gsea_tratadas) <- df_kegg_tratadas_ready$ENTREZID

# 2. Rodar o GSEA para KEGG
set.seed(123)
gsea_kegg_tratadas <- gseKEGG(
  geneList     = vetor_gsea_tratadas,
  organism     = 'hsa',
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose      = FALSE
)

# 3. Traduzir os resultados de volta para o nome normal do gene (facilita leitura)
gsea_kegg_tratadas <- setReadable(gsea_kegg_tratadas, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# =========================================================================
# GRÁFICO DIVIDIDO: UP NO KO vs DOWN NO KO
# =========================================================================

# Alterar a nomenclatura do sinal para ficar autoexplicativo no gráfico
# NES > 0 significa ativado no primeiro termo do contraste (Mock / KO)
gsea_kegg_tratadas@result <- gsea_kegg_tratadas@result %>%
  mutate(.sign = ifelse(NES > 0, "Upregulated in KO", "Downregulated in KO"))

# Gerar o gráfico de pontos (Dotplot) separado
plot_gsea_kegg_dividido <- dotplot(gsea_kegg_tratadas, showCategory = 10, split = ".sign", font.size = 11) +
  facet_grid(~.sign) + 
  labs(
    title = "KEGG Pathways Enrichment (GSEA)",
    subtitle = "HMOX1 KO vs Addback (Treated with Hemin)"
  ) +
  scale_color_gradient(
    low = "#7DA1C4",  # Azul suave para maior significância
    high = "#D67A7A", # Vermelho/rosa suave para menor significância
    name = "p.adjust"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"), 
    strip.background = element_rect(fill = "#f0f0f0", color = "black"), 
    text = element_text(family = "sans")
  )

# Visualizar o gráfico
plot(plot_gsea_kegg_dividido)

# =========================================================================
# GRÁFICO DE DENSIDADE (RIDGEPLOT) - VISUALIZAÇÃO ALTERNATIVA
# =========================================================================
# O ridgeplot mostra a distribuição dos genes dentro de cada via enriquecida
library(ggridges) # Garantir que o pacote está carregado
plot_ridge <- ridgeplot(gsea_kegg_tratadas, showCategory = 15) +
  labs(
    title = "Pathway Expression Distribution",
    subtitle = "Positive density = UP in KO | Negative density = DOWN in KO"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 11, color = "black"),
    text = element_text(family = "sans")
  )

plot(plot_ridge)

# =========================================================================
# 12. ENRICHR VIA R (CONECTANDO À API) E GRÁFICO ESPELHADO
# =========================================================================

# 1. Instalar e carregar pacote do Enrichr
if (!require("enrichR", quietly = TRUE)) install.packages("enrichR")
library(enrichR)
library(dplyr)
library(ggplot2)

# 2. Conectar ao servidor principal (Humano)
setEnrichrSite("Enrichr") 

# DICA: Rodar a linha abaixo abre uma tabela com todos os bancos disponíveis hoje.
# Assim você pode checar se já existe um "KEGG_2026_Human" etc.
dbs <- listEnrichrDbs() 
View(dbs)

# 3. Escolher as bases de dados exatas (Atualize os anos conforme a tabela do dbs)
bancos_desejados <- c("KEGG_2026", "GO_Biological_Process_2026", "Reactome_Pathways_2024", "GO_Molecular_Function_2026", "ChEA_2022")

# 4. Separar as listas de genes usando a sua tabela já tratada (sig_tratadas)
genes_up   <- sig_tratadas %>% filter(logFC > 1) %>% pull(gene)
genes_down <- sig_tratadas %>% filter(logFC < -1) %>% pull(gene)

# 5. Enviar para a API do Enrichr
resultados_up   <- enrichr(genes_up, databases = bancos_desejados)
resultados_down <- enrichr(genes_down, databases = bancos_desejados)

# =========================================================================
# PREPARAÇÃO DA TABELA PARA O GRÁFICO
# =========================================================================

# Função para extrair o Top 10 de um banco específico e juntar UP/DOWN
preparar_tabela_espelhada <- function(enrich_up, enrich_down, nome_banco, top_n = 10) {
  
  # Preparar tabela UP (Barras positivas)
  df_up <- enrich_up[[nome_banco]] %>%
    filter(Adjusted.P.value < 0.05) %>%
    arrange(Adjusted.P.value) %>%
    head(top_n) %>%
    mutate(
      Status = "Upregulated in KO",
      Score = -log10(Adjusted.P.value) # Escore positivo para o gráfico
    )
  
  # Preparar tabela DOWN (Barras negativas)
  df_down <- enrich_down[[nome_banco]] %>%
    filter(Adjusted.P.value < 0.05) %>%
    arrange(Adjusted.P.value) %>%
    head(top_n) %>%
    mutate(
      Status = "Downregulated in KO",
      Score = log10(Adjusted.P.value) # O log de um p-valor < 0.05 já é negativo!
    )
  
  # Unir as tabelas e ordenar para o visual em "escada"
  df_final <- bind_rows(df_up, df_down) %>%
    arrange(Score) %>%
    mutate(Term = factor(Term, levels = unique(Term))) # Trava a ordem pro ggplot não bagunçar
  
  return(df_final)
}

# Escolha qual banco quer plotar agora:
df_plot_reactome <- preparar_tabela_espelhada(resultados_up, resultados_down, "Reactome_Pathways_2024", top_n = 10)
df_plot_kegg <-preparar_tabela_espelhada(resultados_up, resultados_down, "KEGG_2026", top_n = 10)
df_plot_go_bp <-preparar_tabela_espelhada(resultados_up, resultados_down, "GO_Biological_Process_2026", top_n = 10)
df_plot_molecular_function <- preparar_tabela_espelhada(resultados_up, resultados_down, "GO_Molecular_Function_2026", top_n = 20)
df_plot_chea <- preparar_tabela_espelhada(resultados_up, resultados_down, "ChEA_2022", top_n = 10)
# =========================================================================
# GRÁFICO PERSONALIZADO ESTILO LOFI
# =========================================================================

# =========================================================================
# FUNÇÃO PARA PLOTAR ENRICHR (ESTÉTICA LOFI / CLEAN)
# =========================================================================

gerar_grafico_enrichr <- function(df_plot, titulo_banco = "KEGG", subtitulo_exp = "HMOX1 KO vs Addback") {
  
  # Garantir que os termos fiquem na ordem correta do score
  df_plot$Term <- factor(df_plot$Term, levels = df_plot$Term[order(df_plot$Score)])
  
  # ADICIONADO: color = Status dentro do aes()
  grafico <- ggplot(df_plot, aes(x = Score, y = Term, fill = Status, color = Status)) +
    
    # ADICIONADO: Removemos o color = NA e controlamos a espessura da borda com linewidth
    geom_col(width = 0.7, linewidth = 0.6) + 
    
    # Cores de preenchimento (tons pastéis originais)
    scale_fill_manual(values = c(
      "Upregulated in KO"   = "#DCA3A3", # Rosa antigo pastel
      "Downregulated in KO" = "#8FA4B8"  # Azul acinzentado pastel
    )) +
    
    # ADICIONADO: Cores das bordas (tons levemente mais escuros)
    scale_color_manual(values = c(
      "Upregulated in KO"   = "#C48A8A", # Rosa um pouco mais fechado/escuro
      "Downregulated in KO" = "#7188A0"  # Azul um pouco mais fechado/escuro
    )) +
    
    geom_vline(xintercept = 0, color = "#B0B0B0", linewidth = 0.5, linetype = "dashed") +
    
    labs(
      title = paste("Enrichr Pathway Analysis (", titulo_banco, ")", sep = ""),
      subtitle = subtitulo_exp,
      x = "Directional -log10 (Adjusted p-value)",
      y = NULL
    ) +
    
    theme_minimal(base_family = "sans") +
    theme(
      text = element_text(color = "#5c5c5c"), 
      plot.title.position = "plot", 
      plot.title = element_text(face = "bold", size = 15, margin = margin(b = 5), color = "#3d3d3d", hjust = 0),
      plot.subtitle = element_text(size = 11, color = "#8a8a8a", margin = margin(b = 15), hjust = 0),
      axis.text.y = element_text(size = 10.5, color = "#4a4a4a"),
      axis.text.x = element_text(size = 10, color = "#6e6e6e"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      panel.grid.major.x = element_line(color = "#ebebeb", linewidth = 0.5),
      panel.grid.major.y = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
  
  return(grafico)
}

#Gráfico KEGG
plot_kegg <- gerar_grafico_enrichr(df_plot_kegg, titulo_banco = "KEGG 2026")
plot(plot_kegg)

# 2. Gráfico GO Biological Process
plot_go <- gerar_grafico_enrichr(df_plot_go_bp, titulo_banco = "GO Biological Process 2026")
plot(plot_go)

# 3. Gráfico Reactome
plot_reactome <- gerar_grafico_enrichr(df_plot_reactome, titulo_banco = "Reactome 2024")
plot(plot_reactome)

#4. Gráfico Molecular Function
plot_molecular_function <-gerar_grafico_enrichr(df_plot_molecular_function, titulo_banco = "GO Molecular Function 2026")
plot(plot_molecular_function)

#5.Gráfico ChEA
plot_chea <- gerar_grafico_enrichr(df_plot_chea, titulo_banco = "ChEA 2022")
plot(plot_chea)


# =========================================================================
# 13. ENRICHR: CONDIÇÃO BASAL (HMOX1 KO vs ADDBACK - SEM TRATAMENTO)
# =========================================================================

# 1. Extrair as listas baseadas na lógica do contraste (Addback_NT - Mock_NT)
genes_up_ko_basal   <- sig_basal %>% filter(logFC < -1) %>% pull(gene) # Negativo = Maior no KO
genes_down_ko_basal <- sig_basal %>% filter(logFC > 1) %>% pull(gene)  # Positivo = Menor no KO

# 2. Enviar para a API do Enrichr usando os bancos já definidos no passo 12
resultados_up_basal   <- enrichr(genes_up_ko_basal, databases = bancos_desejados)
resultados_down_basal <- enrichr(genes_down_ko_basal, databases = bancos_desejados)

# 3. Preparar as tabelas usando sua função
df_plot_kegg_basal <- preparar_tabela_espelhada(
  resultados_up_basal, 
  resultados_down_basal, 
  "KEGG_2026", 
  top_n = 10
)

df_plot_go_bp_basal <- preparar_tabela_espelhada(
  resultados_up_basal, 
  resultados_down_basal, 
  "GO_Biological_Process_2026", 
  top_n = 10
)

# 4. Gerar os gráficos com a estética personalizada
plot_kegg_basal <- gerar_grafico_enrichr(
  df_plot_kegg_basal, 
  titulo_banco = "KEGG 2026", 
  subtitulo_exp = "Basal Effect: HMOX1 KO vs Addback (Untreated)"
)

plot_go_basal <- gerar_grafico_enrichr(
  df_plot_go_bp_basal, 
  titulo_banco = "GO Biological Process 2026", 
  subtitulo_exp = "Basal Effect: HMOX1 KO vs Addback (Untreated)"
)

# 5. Exibir
plot(plot_kegg_basal)
plot(plot_go_basal)
