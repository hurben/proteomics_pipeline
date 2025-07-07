
import pandas as pd
import sys


def read_file_locations(file_path):

    sample_info_dir = ""
    data_marix_dir = ""
    
    df = pd.read_csv(file_path, sep="\t", header=None, index_col=False)
    r, c = df.shape
  
    for i in range(r):
        if "sample_info" in df.iloc[i, 0]:
            sample_info_dir = df.iloc[i, 1]
        elif "data_matrix" in df.iloc[i, 0]:
            data_marix_dir = df.iloc[i, 1]

    return sample_info_dir, data_marix_dir


if __name__ == "__main__":

    input_file = sys.argv[1]

    sample_info_dir, data_marix_dir = read_file_locations(input_file)

    #manual coding for sample_info_df and data_matrix_df
    #STEP1
    data_matrix_df = pd.read_excel(data_marix_dir, sheet_name=0, skiprows=2)
    data_matrix_df = data_matrix_df.T

    data_matrix_df = data_matrix_df.rename(index={"Assay": "Sample ID"})
    data_matrix_df.rename(index={"Assay": "Sample ID"})
    data_matrix_df = data_matrix_df.reset_index()
    data_matrix_df = data_matrix_df.rename(columns={'index': 'Sample ID'})
    data_matrix_df.columns = data_matrix_df.iloc[0] 
    data_matrix_df = data_matrix_df.iloc[3:]

    #STEP2
    sample_info_df = pd.read_excel(sample_info_dir, sheet_name=0)
    concat_df = sample_info_df.join(data_matrix_df.set_index("Sample ID"), on="Sample ID", how="inner")

    concat_df.to_csv(sys.argv[2], sep="\t", index=False)

