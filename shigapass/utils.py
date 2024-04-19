import pandas as pd

def dict_to_str(kv_dict):
    result = ""
    for key, value in kv_dict.items():
        result += f"{key};{value}\n"
    return result

def write_out(content, outfpath, method="w"):
    with open(outfpath, method) as fout:
        fout.write(content)

def read_file(input_filepath):
    with open(input_filepath, "r") as fin:
        content = fin.read()
    return content

def get_val_from_file(csv_file, sort_column, value_column, delimiter=";"):
    """Get the value from the specified column in the top row that has been sorted via the sort column"""
    try:
        df = pd.read_csv(csv_file, delimiter=delimiter, header=None)
        sorted_df = df.sort_values(by=sort_column)
        top_row_value = sorted_df.iloc[0, value_column]
    except pd.errors.EmptyDataError:
        top_row_value = ""
    return top_row_value
