library(SummarizedExperiment)

message("=========================================")
message("EXTRAÇÃO ROBUSTA DE DADOS (TOP2A)")
message("=========================================")

# 1. Carregar o objeto
fpath <- "data/processed/tcga_brca_prepared_counts.RData"
if (!file.exists(fpath)) stop("Arquivo RData não encontrado!")

message("Lendo arquivo...")
env <- new.env()
nm <- load(fpath, env)[1]
se <- env[[nm]]
message(paste("Objeto carregado:", nm))

# 2. Diagnóstico dos Nomes de Genes
rn <- rownames(se)
message("Exemplo de rownames (primeiros 5):")
print(head(rn))

# Tenta encontrar TOP2A
gene_symbol <- "TOP2A"
gene_id <- NULL

# Estratégia A: O rowname É o símbolo?
if (gene_symbol %in% rn) {
    gene_id <- gene_symbol
    message("✅ Encontrado direto nos rownames!")
}

# Estratégia B: Procurar nos metadados (rowData)
if (is.null(gene_id)) {
    message("Procurando nos metadados (rowData)...")
    rd <- rowData(se)

    # Lista de colunas comuns onde o nome do gene pode estar
    cols_to_check <- c("gene_name", "symbol", "external_gene_name", "gene_symbol", "hgnc_symbol")
    found_col <- NULL

    for (col in cols_to_check) {
        if (col %in% colnames(rd)) {
            message(paste("Checando coluna:", col))
            # Busca exata
            idx <- which(rd[[col]] == gene_symbol)
            if (length(idx) > 0) {
                gene_id <- rownames(se)[idx[1]]
                found_col <- col
                message(paste("✅ Encontrado na coluna:", col, "-> ID:", gene_id))
                break
            }
        }
    }
}

# Estratégia C: Busca "fuzzy" (se for algo como "TOP2A|12345")
if (is.null(gene_id)) {
    message("Tentando busca parcial (grep)...")
    idx <- grep(paste0("^", gene_symbol, "\\|"), rn) # Começa com TOP2A|
    if (length(idx) == 0) idx <- grep(paste0("\\|", gene_symbol, "$"), rn) # Termina com |TOP2A

    if (length(idx) > 0) {
        gene_id <- rn[idx[1]]
        message(paste("✅ Encontrado por busca parcial -> ID:", gene_id))
    }
}

# 3. Extrair e Salvar
if (!is.null(gene_id)) {
    # Tenta pegar assay 'unstranded' ou o primeiro disponivel
    anames <- assayNames(se)
    target_assay <- if ("unstranded" %in% anames) "unstranded" else 1

    counts_vec <- assay(se, target_assay)[gene_id, ]
    clin <- colData(se)

    df_out <- data.frame(
        patient = colnames(se),
        expression = as.numeric(counts_vec),
        days_to_death = if ("days_to_death" %in% names(clin)) clin$days_to_death else NA,
        days_to_last_follow_up = if ("days_to_last_follow_up" %in% names(clin)) clin$days_to_last_follow_up else NA,
        vital_status = if ("vital_status" %in% names(clin)) clin$vital_status else NA
    )

    write.csv(df_out, "survival_data_TOP2A.csv", row.names = FALSE)
    message("\n🎉 SUCESSO! Dados extraídos para 'survival_data_TOP2A.csv'.")
    message("Agora você pode rodar o script Python para gerar o gráfico.")
} else {
    message("\n❌ ERRO FATAL: Gene TOP2A não encontrado de jeito nenhum.")
    message("Colunas disponíveis no rowData:")
    print(colnames(rowData(se)))
    message("Se você ver uma coluna com nomes de genes, edite o script para usar ela.")
}
