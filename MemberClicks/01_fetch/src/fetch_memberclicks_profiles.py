import os
import base64
import requests
import pandas as pd
from dotenv import load_dotenv

def get_api_client_credentials():

    # load variables from .env file
    # should this be in 01_fetch/in ?
    load_dotenv()

    # retrieve client id and client secret
    client_id = os.getenv('CLIENT_ID')
    client_secret = os.getenv('CLIENT_SECRET')

    return client_id,client_secret

def get_access_token(client_id,client_secret):

    # define api client credentials in dictionary
    credentials = f'{
        client_id}:{client_secret
                    }'
    # encode credentials in base64 
    encoded_credentials = base64.b64encode(credentials.encode('utf-8')).decode('utf-8')

    # define url for API request to token endpoint 
    token_url = 'https://cuahsi.memberclicks.net/oauth/v1/token'
    # define headers and payload
    headers = {
        'Authorization': f'Basic {encoded_credentials}',
        'Content-Type': 'application/x-www-form-urlencoded'
    }
    payload = {
        'grant_type': 'client_credentials'
    }
    try:
        # need to send POST request to token API endpoint with encoded data
        response = requests.post(token_url,headers=headers,data=payload)
        response.raise_for_status()

        if response.status_code == 200:
            # extract the json response
            json_response = response.json()
            # retrieve the access token if JSON response is not empty
            if json_response is not None:
                return json_response.get('access_token')
            else:
                print('Error in obtaining access token.')
                return None
    # error handling
    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except requests.exceptions.ConnectionError as conn_err:
        print(f"Connection error occurred: {conn_err}")
    except requests.exceptions.Timeout as timeout_err:
        print(f"Timeout error occurred: {timeout_err}")
    except requests.exceptions.RequestException as req_err:
        print(f"An error occurred: {req_err}")

def get_memberclicks_profiles(access_token):

    # define headers for API request
    headers = {
        'Authorization': f'Bearer {access_token}'
        }
    try:
        # initialize list of profiles (contained in dictionaries)
        profiles = []

        # endpoint for first page of profiles 
        page_profile_url = 'https://cuahsi.memberclicks.net/api/v1/profile?pageSize=100&pageNumber=1'

        # create loop to build list of profiles
        while page_profile_url:
            # need to send GET request to profile API endpoint
            response = requests.get(page_profile_url, headers=headers)
            response.raise_for_status()
            if response.status_code == 200:
                # extract the json response
                json_response = response.json()
                # retrieve the profiles if JSON response is not empty
                if json_response is not None:
                    # build (add) list of profiles
                    profiles += json_response.get('profiles')
                    # check if there is a next page
                    if json_response.get('nextPageUrl') is not None:
                        page_profile_url = json_response.get('nextPageUrl')
                    else:
                        # this will break the while loop
                        page_profile_url = None
                else:
                    # break loop if no profiles are found
                    print('Error in obtaining profiles.')
                    break
        # convert list of dictionaries to dataframe and return
        return pd.DataFrame(profiles)
            
    # error handling
    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except requests.exceptions.ConnectionError as conn_err:
        print(f"Connection error occurred: {conn_err}")
    except requests.exceptions.Timeout as timeout_err:
        print(f"Timeout error occurred: {timeout_err}")
    except requests.exceptions.RequestException as req_err:
        print(f"An error occurred: {req_err}")

def save_output_file(df,out_file):
    
    df.to_csv(out_file,index=False)

def main(out_file):
    print('Loading...')
    # retrieve client id and client secret
    client_id,client_secret = get_api_client_credentials()
    print(f'Client id = {client_id}')
    # retrieve access token
    access_token = get_access_token(client_id,client_secret)
    print('Access token retrieved.')
    # get profiles
    df_profiles = get_memberclicks_profiles(access_token)

    save_output_file(df_profiles,out_file)

if __name__ == '__main__':

    # inputs from snakefile
    out_filename = snakemake.output['out_filename']

    # main function
    main(out_filename)


