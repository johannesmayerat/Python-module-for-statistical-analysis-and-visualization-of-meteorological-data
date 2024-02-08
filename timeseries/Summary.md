df.describe()             print statistics
df.loc['...']                   select row  by df index, also list possible ['...','....']
df.iloc[-1,1:5]             select rows and cols by numeric index
df.index                     return dataframe index ( set by set_index)
df.idxmax()               return df index of highest value in each column

columns = [col for col in df.columns if col[:2] == 'AT']

df.isna().sum()
df = df.dropna(subset=['A','B'], axis=0)
df.replace('?',np.NaN)

df['price'].astype('int')
df['Date'] = pd.to_datetime(df['Date']).set_index('Date')    # df mit .copy() erstellt

vallist = df_aut.index[df_aut['solar_generation_actual'] == 37.0].tolist()
df_aut.loc[vallist]

df['...'].value_counts()
df['...'].unique()

df_monthly = df.reset_index().resample('M',on='Date').mean()
# nicht mit df index möglich, daher voher reset_index

df_monthly.rolling(window=12, center=True).mean()

# for seasonal boxplot; if df index is timestamp obj (and not simply a object)
df_aut_monthly['month'] = df_aut_monthly.index.month

df['...'].plot(kind='line')               kind = line, hist, area, bar, barh, pie, box  

df.sort_value(['...'], ascending=True, axis=0, inplace=True)
