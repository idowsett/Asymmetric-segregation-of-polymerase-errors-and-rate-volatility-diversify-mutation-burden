

# mutations = [
# 248	,200,	255	,328,	288,	279,
# 285 ,231,	256	,251,	262,	301,
# 255	,287,	298	,270,	277,	316	,308,	356,	307,	306,
# 270	,289,	250	,303,	155,	334,
# 231	,236,	298	,302,	283,
# 117	,323,	312	,288,	295,	335,
# 252	,282,	319	,296,	284,	254,	224,	253,	264,	280,	232,	219
# ]
# print(len(division))
# print(len(mutations))
# dictionary = {
# "L_151": 2
# ,"L_151" : 3
# ,"L_151" : 4
# ,"L_151" : 6
# ,"L_151" : 7
# ,"L_151" : 8
# ,"L_153" : 10
# ,"L_153" : 11
# ,"L_153" : 12
# ,"L_153" : 13
# ,"L_153" : 14
# ,"L_153" : 15
# ,"L_156" : 2
# ,"L_156" : 3
# ,"L_156" : 4
# ,"L_156" : 5
# ,"L_156" : 6
# ,"L_156" : 7
# ,"L_156" : 8
# ,"L_156" : 10
# ,"L_156" : 11
# ,"L_156" : 12
# ,"L_157" : 2
# ,"L_157" : 3
# ,"L_157" : 4
# ,"L_157" : 5
# ,"L_157" : 6
# ,"L_157" : 7
# ,"L_158" : 2
# ,"L_158" : 3
# ,"L_158" : 5
# ,"L_158" : 6
# ,"L_158" : 8
# ,"L_160" : 2
# ,"L_160" : 3
# ,"L_160" : 4
# ,"L_160" : 6
# ,"L_160" : 8
# ,"L_160" : 10
# ,"L_162" : 2
# ,"L_162" : 3
# ,"L_162" : 4
# ,"L_162" : 5
# ,"L_162" : 7
# ,"L_162" : 8
# ,"L_162" : 11
# ,"L_162" : 12
# ,"L_162" : 13
# ,"L_162" : 14
# ,"L_162" : 15
# ,"L_162" : 16}


data = {'division':[2,3,4,6,7,8,10,11,12,13,14,15,2,3,4,5,6,7,8,10,11,12,2,3,
4,5,6,7,2,3,5,6,8,3,4,6,8,10,2,3,4,5,7,8,11,12,13,14,15,16],
# 'mutations':[248,200,255,328,288,279,285,231,256,251,262,301,255,287,298,270,277,316,308,356,307,306,270,289,250,303,
# 155,334,231,236,298,302,283,117,323,312,288,295,335,252,282,319,296,284,254,224,253,264,280,232,219]}
'mutations' : [
248	,200,	255	,328,	288,	278,
284 ,228,	256	,251,	262,	301,
253	,288,	297	,270,	278,	315	,306,	356,	307,	306,
269	,288,	250	,303,	155,	334,
229	,235,	298	,302,	283,
323,	312	,288,	288,	332,
252	,283,	319	,298,	281,	254,	224,	252,	263,	280,	232,	219
]}

l151x = [2, 3, 4, 6, 7, 8]
l151y = [248,200,255,328,288,278]

l153x = [10, 11, 12, 13, 14, 15]
l153y = [284,228,256,251,262,301]

l156x = [2, 3, 4, 5, 6, 7, 8, 10, 11, 12]
l156y = [253,288,297,270,278,315,306,356,307,306]


l157x = [2, 3, 4, 5, 6, 7]
l157y = [269,288,250,303,155,334]


l158x = [2, 3, 5, 6, 8]
l158y = [229,235,298,302,283]


l160x = [3, 4, 6, 8, 10]
l160y = [323,312,288,288,332]


l162x = [2, 3, 4, 5, 7, 8, 11, 12, 13, 14, 15, 16]
l162y = [252,283,319,298,281,254,224,252,263,280,232,219]




import pandas as pd
df = pd.DataFrame(data, columns = ['division','mutations'])

from scipy import stats
import plotly
import plotly.graph_objs as go
import plotly.figure_factory as ff
spearman_coef, ps_value = stats.spearmanr(df["division"], df["mutations"])
print("Spearman Correlation Coefficient: ", spearman_coef, "and a P-value of:", ps_value)

# Plotly mapping of relationship, all lineages seperately
# Title = "Correlation Between Age of Mother and Number of Mutations <br> Pearson Correlation Coefficient: " + str(round(pearson_coef,5)) + "\n, p = " + str(round(p_value, 3))
Title = "All lineages"
lin151 = go.Scatter(
    x = l151x,    y = l151y,    mode = 'lines',    line = dict(width = 6), marker=dict(        size=20,        symbol='circle',        color= 'black',        opacity=1 ))
lin153 = go.Scatter(
    x = l153x,    y = l153y,    mode = 'lines',    line = dict(width = 6), marker=dict(        size=8,        symbol='circle',        color= 'hotpink',        opacity=1 ))
lin156 = go.Scatter(
    x = l156x,    y = l156y,    mode = 'lines',    line = dict(width = 6), marker=dict(        size=8,        symbol='circle',        color= 'green',        opacity=1 ))
lin157 = go.Scatter(
    x = l157x,    y = l157y,    mode = 'lines',    line = dict(width = 6), marker=dict(        size=8,        symbol='circle',        color= 'blue',        opacity=1 ))
lin158 = go.Scatter(
    x = l158x,    y = l158y,    mode = 'lines',    line = dict(width = 6), marker=dict(        size=8,        symbol='circle',        color= 'orange',        opacity=1 ))
lin160 = go.Scatter(
    x = l160x,    y = l160y,    mode = 'lines',    line = dict(width = 6), marker=dict(        size=8,        symbol='circle',        color= 'purple',        opacity=1 ))
lin162 = go.Scatter(
    x = l162x,    y = l162y,    mode = 'lines',    line = dict(width = 6), marker=dict(        size=8,        symbol='circle',        color= 'red',        opacity=1 ))

data1 = [
lin151,
lin153,
lin156,
lin157,
lin158,
lin160,
lin162
]
layout = go.Layout(
    title=Title,
    font=dict(family='Arial, sans-serif', size=10, color='black'),
    yaxis=dict(title='Mutations', autorange=False, range=[0, 400], titlefont=dict(family='Arial, sans-serif', size=14, color='black')),
    xaxis=dict(title='Division Number', autorange=False, range=[1, 17], titlefont=dict(family='Arial, sans-serif', size=14, color='black')),
    autosize=False,
    width=700,
    height=700,
    margin=dict(l=50, r=10, b=50, t=90, pad=4), paper_bgcolor='#ffffff', plot_bgcolor='#ffffff'
    )
# layout = go.Layout(title = Title)
fig = go.Figure(data=data1, layout=layout)
plotly.offline.plot(fig, auto_open=True)
