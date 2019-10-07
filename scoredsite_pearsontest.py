import plotly
import plotly.graph_objs as go
import plotly.figure_factory as ff
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats


data = {'scoredsites':[10541044,10541044,10541044,10541044,10541044,10541044,
10636360,10636360,10636360,10636360,10636360,10636360,
9521963,9521963,9521963,9521963,9521963,9521963,9521963,9521963,9521963,9521963,
10594546,10594546,10594546,10594546,10594546,10594546,
9818623,9818623,9818623,9818623,9818623,
9679326,9679326,9679326,9679326,9679326,
8578576,8578576,8578576,8578576,8578576,8578576,8578576,8578576,8578576,8578576,8578576,8578576],
'mutations':
[248	,200,	255	,328,	288,	278,
284 ,228,	256	,251,	262,	301,
253	,288,	297	,270,	278,	315,   306,	356,	307,	306,
269	,288,	250	,303,	155,	334,
229	,235,	298	,302,	283,
323 ,312,   288 ,288,	332,
252	,283,	319	,298,	281,	254,	224,	252,	263,	280,	232,	219]}



# Convert data to dataframe
df = pd.DataFrame(data, columns = ['scoredsites','mutations'])
pearson_coef, p_value = stats.pearsonr(df["scoredsites"], df["mutations"])
print("Pearson Correlation Coefficient: ", pearson_coef, "and a P-value of:", p_value)

Title = "Number of Scored Sites and Mutations Correlation"
stats = "r = " + str(round(pearson_coef,5)) + "\np = " + str(round(p_value, 3))

# Generating scatterplot:
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
sns.set(style="white", palette="bright", color_codes=True, rc={"lines.linewidth": 1.0})
ax = sns.scatterplot(x='scoredsites', y='mutations', data=df)
plt.xlabel('Scored Sites ( * 1e7)',fontsize=11)
plt.ylabel('Mutations', fontsize=11)
ax.get_xaxis().get_offset_text().set_position((0.1,0))
plt.title(Title
# , loc='left'
)
plt.text(8528576, 150, stats, fontsize=12)
plt.savefig(f"ScoredSiteCorrelation.pdf",transparent = False,bbox_inches='tight')
plt.show()
