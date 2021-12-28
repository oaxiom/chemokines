
from ete3 import Tree, TreeStyle

t = Tree('single_representative_envs.dnd')

ts.show_leaf_name = True
ts.mode = "c"

t.render('envelopes.pdf', tree_style=ts)
