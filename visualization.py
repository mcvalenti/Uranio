import pandas as pd
import plotly.graph_objects as go


columns = ['x','y','z','vx','vy','vz','m']
df1 = pd.read_csv('output_Files/orbit1.csv', sep='\t',  names=columns, header=None, index_col=False)
df1['orbit']=1

df2 = pd.read_csv('output_Files/orbit2.csv', sep='\t',  names=columns, header=None, index_col=False)
df2['orbit']=2

df3 = pd.read_csv('output_Files/orbit3.csv', sep='\t',  names=columns, header=None, index_col=False)
df3['orbit']=3

df4 = pd.read_csv('output_Files/orbit4.csv', sep='\t',  names=columns, header=None, index_col=False)
df4['orbit']=4

df_total = pd.concat([df1, df2, df3, df4], ignore_index=True)
print(df_total.info())


# PLOT
fig = go.Figure()
for nombre, grupo in df_total.groupby('orbit'):
    fig.add_trace(go.Scatter(
        x=grupo['x'],
        y=grupo['y'],
        mode='lines+markers',
        name=nombre  # nombre del archivo como etiqueta
    ))

fig.update_layout(
    title='Trajectory [km]',
    xaxis_title='X [km]',
    yaxis_title='Y [km]',
    width=500,
    height=500,
)
fig.show()

