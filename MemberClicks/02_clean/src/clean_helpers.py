import pandas as pd

def rename_columns(df,new_columns):
    # this function is not used at the moment 

    # rename dataframe columns in vectorized way
    df.rename(columns=dict(zip(df.columns, new_columns)), inplace=True)

    return df

def remove_columns(df,columns_to_remove):

    df.drop(columns=columns_to_remove,inplace=True)

    return df

def remove_brackets_from_columns(df):
    # some MC exports retain square brackets around column headers
    for col in df.columns:
        new_col = col.replace('[', '').replace(']', '')
        df.rename(columns={col:new_col},inplace=True)

    return df

def modify_institution_names(df,institution_field):

    # remove special unicode characters (e.g. U of Hawaii)
    df[institution_field] = df[institution_field].str.encode('utf-8').str.decode('ascii', 'ignore')

    # create dictionary of instituion values to be updated
    # note we are treating UNC system as one institution (as opposed to segmenting by campus (e.g. Chapel Hill))
    replace_vals = {
        'University of North Carolina': 'University of North Carolina System'
        }
    
    # update institution names
    df[institution_field] = df[institution_field].replace(replace_vals)

    return df