import pandas as pd
import clean_helpers

def read_input_file(in_file):
    
    df = pd.read_csv(in_file)

    return df

def add_payment_line(df):

    # group by invoice ID and assign a sequential number to each payment line
    df['Payment Line'] = df.groupby('Invoice ID').cumcount() + 1
    
    return df

def save_output_file(df,out_file):
    
    df.to_csv(out_file,index=False)
    
def main(in_file,out_file):

    columns_to_remove = ['Transaction Type',
                         'Payment Method',
                         'Cardholder Name',
                         'Transaction Status',
                         'Response Text',
                         'Reference Number',
                         'Revenue Account',
                         'Credit Account',
                         'Debit Account',
                         'First Name',
                         'Last Name']
  
    df = read_input_file(in_file) 

    print(df.columns)

    df_cleaned = clean_helpers.remove_columns(df,columns_to_remove)

    df_cleaned = add_payment_line(df_cleaned)

    df_cleaned = clean_helpers.modify_institution_names(df_cleaned,'Company')

    save_output_file(df_cleaned,out_file)

if __name__ == '__main__':

    # inputs from snakefile
    in_filename = snakemake.input['in_filename']
    out_filename = snakemake.output['out_filename']

    # main function
    main(in_filename,out_filename)