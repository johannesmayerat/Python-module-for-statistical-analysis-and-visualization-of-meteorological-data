df.describe()             print statistics
df.loc['...']                   select row  by df index, also list possible ['...','....']
df.iloc[-1,1:5]             select rows and cols by numeric index
df.index                     return dataframe index ( set by set_index)
df.idxmax()               return df index of highest value in each column
df['...'].value_counts()
df['...'].unique()

# select specific columns
columns = [col for col in df.columns if col[:2] == 'AT']

# remove nans
df.isna().sum()
df = df.dropna(subset=['A','B'], axis=0)
df.replace('?',np.NaN)

# drop columns
df.drop(['B', 'C'], axis=1)

# drop rows by index
df.drop([0, 1, 7])

# convert column type
df['price'].astype('int')
df['Date'] = pd.to_datetime(df['Date']).set_index('Date')    # df mit .copy() erstellt

# Sort dataframe
df.sort_value(['...'], ascending=True, axis=0, inplace=True)

# index search
vallist = df_aut.index[df_aut['solar_generation_actual'] == 37.0].tolist()
df_aut.loc[vallist]

# wenn Date nicht df index ist
df_monthly = df.resample('M',on='Date').mean()
# wenn Date index ist
df_aut_monthly = df_aut.resample('M').mean()

# apply moving avaergae
df_monthly.rolling(window=12, center=True).mean()

# for seasonal boxplot; if df index is timestamp obj (and not simply a object)
df_aut_monthly['month'] = df_aut_monthly.index.month

# Plot
df['...'].plot(kind='line')               kind = line, hist, area, bar, barh, pie, box  
