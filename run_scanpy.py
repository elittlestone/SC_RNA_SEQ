import scanpy as sc
import anndata as ad
import argparse

def main():
    parser = argparse.ArgumentParser(description = "Scanpy scRNA-seq processing script")
    parser.add_argument("--sample_id", type = str, required = True,
                        help = "Name of sample")
    parser.add_argument("--feature_matrix_file", type = str, required = True,
                        help = "Path to filtered_feature_bc_matrix file")
    parser.add_argument("--processed_h5", type = str, required = True,
                        help = "Processed h5 file")
    args = parser.parse_args()
    matrix_file = args.feature_matrix_file
    outfile = args.processed_h5
    sample_id = args.sample_id

def read_10x_file_into_adata(matrix_file, sample_id):
    adatas = {}
    sample_adata = sc.read_10x_h5(matrix_file)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

    adata = ad.concat(adatas, label = "sample")
    adata.obs_names_make_unique()
    print(adata.obs["sample"].value_counts())
    adata

if __name__ == "__main__":
    main()
