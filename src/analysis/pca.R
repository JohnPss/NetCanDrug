# Instalar pacotes se necessario
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("edgeR", quietly = TRUE)) {
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("edgeR")
}
if (!require("SummarizedExperiment", quietly = TRUE)) {
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("SummarizedExperiment")
}

library(ggplot2)
library(edgeR)
library(SummarizedExperiment)

# 1. Configurar Caminhos
input_file <- "data/processed/tcga_brca_prepared_counts.RData"
output_dir <- "results/figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("=========================================")
message("GERANDO FIGURE S1: PCA (R SCRIPT FIXED)")
message("=========================================")

# 2. Carregar Dados
message(paste("Lendo:", input_file))
env <- new.env()
nm <- load(input_file, env)[1]
se_object <- env[[nm]]

message(paste("Objeto carregado:", nm))
message(paste("Classe do objeto:", class(se_object)[1]))

# --- CORREÇÃO DO ERRO ---
# Extrair a matriz de contagens do objeto SummarizedExperiment
# Tenta encontrar o nome correto do assay (counts, unstranded, etc)
assay_names <- assayNames(se_object)
message("Assays disponíveis no objeto: ", paste(assay_names, collapse = ", "))

# Pega o primeiro assay disponível (geralmente é o de contagem bruta)
# Se houver 'unstranded' (comum no GDC mais novo), prefere ele.
if ("unstranded" %in% assay_names) {
    counts_matrix <- assay(se_object, "unstranded")
    message("Usando assay: 'unstranded'")
} else if ("counts" %in% assay_names) {
    counts_matrix <- assay(se_object, "counts")
    message("Usando assay: 'counts'")
} else {
    counts_matrix <- assay(se_object, 1)
    message(paste("Usando o primeiro assay disponível:", assay_names[1]))
}

message(paste("Dimensões da matriz de contagem:", paste(dim(counts_matrix), collapse = " x ")))

# 3. Preparar Dados para PCA
# Criar objeto DGEList (formato nativo do edgeR)
dge <- DGEList(counts = counts_matrix)

# Filtrar genes com baixa expressão
# Mantém genes com CPM > 1 em pelo menos 10 amostras (ajuste conforme necessário)
keep <- rowSums(cpm(dge) > 1) >= 10
dge_filt <- dge[keep, , keep.lib.sizes = FALSE]
message(paste("Genes retidos após filtro:", nrow(dge_filt)))

# Normalizar (TMM + logCPM)
dge_filt <- calcNormFactors(dge_filt)
logcpm <- cpm(dge_filt, log = TRUE)

# 4. Calcular PCA
message("Calculando Componentes Principais (PCA)...")
# Transpor (amostras nas linhas) e escalar
pca_res <- prcomp(t(logcpm), scale. = TRUE)

# Variância explicada
var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

# 5. Identificar Grupos pelo Barcode TCGA
samples <- rownames(pca_res$x)
# Barcode formato: TCGA-XX-XXXX-01A...
# O 4º campo (01, 11) indica o tipo de tecido
code <- substr(unlist(lapply(strsplit(samples, "-"), function(x) x[4])), 1, 2)

group <- ifelse(as.numeric(code) < 10, "Tumor",
    ifelse(as.numeric(code) >= 10, "Normal", "Outro")
)

df_plot <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Group = group
)

# Manter apenas Tumor e Normal
df_plot <- df_plot[df_plot$Group %in% c("Tumor", "Normal"), ]
message(paste("Amostras finais no plot:", nrow(df_plot)))

# 6. Gerar Figura
message("Gerando gráfico...")

p <- ggplot(df_plot, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(alpha = 0.6, size = 2) +
    stat_ellipse(level = 0.95, show.legend = FALSE, type = "norm", linetype = 2) +
    scale_color_manual(values = c("Normal" = "#377EB8", "Tumor" = "#E41A1C")) +
    theme_bw() +
    labs(
        title = "Figure S1: PCA of Gene Expression (TCGA-BRCA)",
        subtitle = paste("Total samples:", nrow(df_plot)),
        x = paste0("PC1 (", var_explained[1], "% variance)"),
        y = paste0("PC2 (", var_explained[2], "% variance)"),
        color = "Sample Type"
    ) +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "top"
    )

outfile <- file.path(output_dir, "figure_S1_pca.png")
ggsave(outfile, plot = p, width = 8, height = 6, dpi = 300)

message(paste("✅ Sucesso! Imagem salva em:", outfile))
