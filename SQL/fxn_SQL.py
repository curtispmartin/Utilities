# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 10:42:05 2021

@author: MartinC
"""

### utilities
# import sys

### data 
import numpy as np
import pandas as pd

### IO
import pyodbc
import pickle


##### CREATE FUNCTION FOR PULLING DATA FROM SQL SERVER
#----------------------------------------------------------------------------#
### define function... table arg optional
def pull_SQL(db, server, driver='SQL Server', get=None, cols=None, data_type='table', check_nrows=True):
    '''
    This function simplifies the process of pulling data from SQL servers...
    
    db: Name of database (string)
    server: Name of server (string)
    driver: Input for connecting to db; don't change unless absolutely necessary (string)
    get: Name of table or view (string); if None, function will ask you upon pulling table names
    cols: Columns you want to pull (list); requires a priori knowledge 
    data_type: 'table' or 'view' (string); generally want the former
    check_nrows: Checks how many rows are in table (Boolean); if above a certain amount, asks user how many to pull (in order to save RAM if exploratory); set to False if you need complete table
    '''

### create db connection
    conn = pyodbc.connect(
        f'Driver={driver};'
        f'Server={server};'
        f'Database={db};'
        'Trusted_Connection=yes;' # will need more work for untrusted connections
        ) 

### if user wants a table...
    if data_type == 'table':    
        
### display table names for user selection
        sql_TYPE = '''
            SELECT TABLE_NAME
            FROM [{}].INFORMATION_SCHEMA.TABLES
            WHERE TABLE_TYPE = 'BASE TABLE'
            '''.format(db)
        col_name = 'TABLE_NAME' # for pulling right data later...

### if user wants a view...
    elif data_type == 'view':   

### display names for user selection
        sql_TYPE = 'SELECT * FROM sys.views'
        col_name = 'name' # for pulling right data later...

### get names of given data types (table, view, etc.) in database        
    df_TYPE = pd.read_sql(con=conn, sql=sql_TYPE)
    
### return error message & break out of frame if empty... NEED A BETTER WAY TO CHECK CONNECTION ISSUES
    if df_TYPE.empty:
        print('\nNo data found -- Check server connection or permissions. Exiting...')
        conn.close()
        return

### if table name provided, check that it's in database
    if get:
        if get not in list(df_TYPE[col_name]):
            print('\nRequested table not in database. Exiting...')
            conn.close()
            return
        
### otherwise pull table interactively
    else:
        get = df_TYPE.loc[int(input(f'{df_TYPE}\n\nSelect Index: '))][col_name]

### pull columns of interest if indicated
    if cols:
        col_format = ', '.join([f'[{col}]' if ' ' in col else col for col in cols]) # in case any column names have spaces
        
    else:
        col_format = '*'

### sql query        
    sql = '''
        SELECT {} FROM [{}].[dbo].[{}]
        '''.format(col_format, db, get)

### IC1: check number of rows in table (can turn off using `check_nrows=False`)
    if check_nrows:
        
### get length of data... if large, warn user
        sql_NROWS ='''
            SELECT COUNT(*) FROM [{}].[dbo].[{}]
            '''.format(db, get)
        nrows = pd.read_sql(con=conn, sql=sql_NROWS).loc[0].squeeze()
    
### set limit for "big" data (nrows)
        lim_nrows = 100000
        
### if large, ask user how to proceed
        if nrows > lim_nrows:

### constrain number of attempts
            n_tries = 0

### IC2: create loop for allowing multiple tries 
            while True:
                IC1 = input(f'{nrows} rows... Do you want all of it (y/n)? ')
            
### if user wants to load whole thing, g'head...        
                if IC1 == 'y':
                    sql = '''
                        SELECT {} FROM [{}].[dbo].[{}]
                        '''.format(col_format, db, get)
                    break
    
### if not, does user want to load subset?
                elif IC1 == 'n':
                    req_nrows = int(input('How many rows would you like to pull? '))
                    sql = '''
                        SELECT TOP {} {} FROM [{}].[dbo].[{}]
                        '''.format(req_nrows, col_format, db, get)
                    break

### if answer unknown, try again (three times total)    
                else:
                    if n_tries < 2:
                        print('Unknown response... Please try again.')
                        n_tries += 1
                        continue
                    else:
                        print('Unknown response... Exiting.')
                        conn.close()
                        return

### pull requested data
    df = pd.read_sql(con=conn, sql=sql)
    print(f'\n{df.head()}')

### close db connection
    conn.close()

### return requested data only
    return(df)
