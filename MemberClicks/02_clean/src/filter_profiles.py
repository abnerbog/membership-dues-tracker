import pandas as pd
import clean_helpers

def read_input_file(in_file):
    
    df = pd.read_csv(in_file)

    return df

def remove_test_profiles(df):

    # first remove test profiles by filtering rows containing 'test' in the [Organization] field
    df_real = df[~df['[Organization]'].str.contains('test', case=False, na=False)]

    return df_real

def get_institution_profiles(df):

    # note we do not include insititutions with a member type of 'Inactive Members'
    institutions_member_type = [
        'University Member',
        'Primarily Undergraduate Institution',
        'Non-profit  Affiliate Members',
        'International Affiliate Members'
    ] 

    # define list of profile attributes to keep; we can rename columns elsewhere
    institution_fields = ['[Organization]','[Member Type]','[Last Modified Date]','[Profile ID]']

    # filter representative profiles
    df_institutions = df[df['[Member Type]'].isin(institutions_member_type)][institution_fields]

    return df_institutions

def get_representative_profiles(df):

    # note we do not include representatives with a member type of 'Inactive Representatives'
    representatives_member_type = [
        'University Representative',
        'Primarily Undergraduate Institution Representative',
        'Non-profit  Affiliate Representative',
        'International Affiliate Representative'
    ]

    # define list of profile attributes to keep; we can rename columns elsewhere
    representative_fields = ['[Name | Last]','[Name | First]','[Organization]','[Email | Primary]','[Member Type]','[Last Modified Date]']

    # filter representative profiles
    df_representatives = df[df['[Member Type]'].isin(representatives_member_type)][representative_fields]

    return df_representatives

def save_output_file(df,out_file):
    
    df.to_csv(out_file,index=False)

def main(in_file,out_file_representatives,out_file_institutions):

    # read input file
    df = read_input_file(in_file)

    # remove test profiles
    df_real = remove_test_profiles(df)

    # modify institution name
    df_real = clean_helpers.modify_institution_names(df_real,'[Organization]')

    # get representative profiles
    df_representatives = get_representative_profiles(df_real)
    df_representatives = clean_helpers.remove_brackets_from_columns(df_representatives)

    # get institution profiles 
    df_institutions = get_institution_profiles(df_real)
    df_institutions = clean_helpers.remove_brackets_from_columns(df_institutions)

    # save profiles to file 
    save_output_file(df_representatives,out_file_representatives)
    save_output_file(df_institutions,out_file_institutions)

if __name__ == '__main__':

    # inputs from snakefile
    in_filename = snakemake.input['in_filename']
    out_filename_representatives = snakemake.output['out_filename_representatives']
    out_filename_institutions = snakemake.output['out_filename_institutions']

    # main function
    main(in_filename,out_filename_representatives,out_filename_institutions)    