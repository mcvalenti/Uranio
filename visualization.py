import numpy as np
import pandas as pd
import plotly.graph_objects as go

GM=3.98600448e14 # [m*3/s*2]

columns = ['x','y','z','vx','vy','vz','m']

df0 = pd.read_csv('output_Files/orbit0.csv', sep='\t',  names=columns, header=None, index_col=False)
df0['orbit']=0

df1 = pd.read_csv('output_Files/orbit1.csv', sep='\t',  names=columns, header=None, index_col=False)
df1['orbit']=1

df2 = pd.read_csv('output_Files/orbit2.csv', sep='\t',  names=columns, header=None, index_col=False)
df2['orbit']=2

df3 = pd.read_csv('output_Files/orbit_3.csv', sep='\t',  names=columns, header=None, index_col=False)
df3['orbit']=3

df4 = pd.read_csv('output_Files/orbit_4.csv', sep='\t',  names=columns, header=None, index_col=False)
df4['orbit']=4

df_total = pd.concat([df0, df1, df2], ignore_index=True) #  df1, df2, df3, df4
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
    title='Trajectory [m]',
    xaxis_title='X [m]',
    yaxis_title='Y [m]',
    width=500,
    height=500,
)
fig.show()

# df0["r"] = np.sqrt(df0.x**2 + df0.y**2 + df0.z**2)  # Position
# df0["v2"] = df0.vx**2 + df0.vy**2 + df0.vz**2       # Velocity
# df0["energy"] = 0.5*df0.v2 - GM/df0.r              # Energy [J/kg]
# df0["energy_diff"] = df0["energy"] - df0["energy"].mean()

# # PLOT
# fig1 = go.Figure()
# fig1.add_trace(go.Scatter(
#         x=df0.x,
#         y=df0.energy_diff,
#         mode='lines+markers',
#         name=nombre  # nombre del archivo como etiqueta
#     ))

# fig1.update_layout(
#     title='Energy [J/kg]',
#     xaxis_title='X [m]',
#     yaxis_title='Energy [J/kg]'
# )
# fig1.show()
