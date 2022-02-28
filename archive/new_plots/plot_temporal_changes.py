fig_snp, ax_snp = plt.subplots(figsize=(10,6))

colormap = cmx.get_cmap('jet', 6)
colors = [colormap(x) for x in numpy.array([x for x in range(0,6)])/6.0]

# Unrelated mother/infant @
# Unrelated adult @
# Adult-adult @
# Mother-mother @
# Infant-infant @
# Mother-infant @

# Plot SNP change distribution

ax_snp.set_xscale('log')
ax_snp.set_yscale('log')
ax_snp.set_ylabel('Fraction comparisons $\geq n$', fontsize=11)
ax_snp.set_xlabel('# SNP changes', fontsize=11)

ax_snp.spines['top'].set_visible(False)
ax_snp.spines['right'].set_visible(False)
ax_snp.get_xaxis().tick_bottom()
ax_snp.get_yaxis().tick_left()

color_i = 0

# Within-host, adult
counts = []
for tp_pair in pooled_snp_change_distribution['hmp'].keys():
	counts += pooled_snp_change_distribution['hmp'][tp_pair]

xs, ns = calculate_unnormalized_survival_from_vector(counts)
mlabel = 'HMP: Adult 6 months' + (' (n=%d)' % ns[0])
ax_snp.step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label=mlabel, where='pre',zorder=4)
ymin = 1.0/ns[0]
ymax = 1.3
ax_snp.set_ylim([ymin,ymax])
color_i += 1

# Unrelated adults
counts = []
for tp_pair in pooled_between_snp_change_distribution['hmp'].keys():
	counts += pooled_between_snp_change_distribution['hmp'][tp_pair]

xs, ns = calculate_unnormalized_survival_from_vector(counts)
ax_snp.step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label="Unrelated adults", where='pre',zorder=4)
color_i += 1

# Within-host, infant-infant
counts = []
for cohort in mi_cohorts:
	for tp_pair in pooled_snp_change_distribution[cohort].keys():
		tp1, tp2 = tp_pair
		if tp1[0] == 'I' and tp2[0] == 'I':
			counts += pooled_snp_change_distribution[cohort][tp_pair]

xs, ns = calculate_unnormalized_survival_from_vector(counts)
mlabel = "Infant-infant" + (' (n=%d)' % ns[0])
ax_snp.step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label=mlabel, where='pre',zorder=4)
color_i += 1

# Within-host, mother-infant
counts = []
for cohort in mi_cohorts:
	for tp_pair in pooled_snp_change_distribution[cohort].keys():
		tp1, tp2 = tp_pair
		if (tp1[0] == 'I' and tp2[0] == 'M') or (tp1[0] == 'M' and tp2[0] == 'I'):
			counts += pooled_snp_change_distribution[cohort][tp_pair]

xs, ns = calculate_unnormalized_survival_from_vector(counts)
mlabel = "Mother-infant" + (' (n=%d)' % ns[0])
ax_snp.step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label=mlabel, where='pre',zorder=4)
color_i += 1

# Within-host, mother-mother
counts = []
for cohort in mi_cohorts:
	for tp_pair in pooled_snp_change_distribution[cohort].keys():
		tp1, tp2 = tp_pair
		if tp1[0] == 'M' and tp2[0] == 'M':
			counts += pooled_snp_change_distribution[cohort][tp_pair]

xs, ns = calculate_unnormalized_survival_from_vector(counts)
mlabel = "Mother-mother" + (' (n=%d)' % ns[0])
ax_snp.step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label=mlabel, where='pre',zorder=4)
color_i += 1

# Unrelated mother/infant
counts = []
for cohort in mi_cohorts:
	for tp_pair in pooled_between_snp_change_distribution[cohort].keys():
		counts += pooled_between_snp_change_distribution[cohort][tp_pair]

xs, ns = calculate_unnormalized_survival_from_vector(counts)
ax_snp.step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label="Unrelated mothers/infants", where='pre',zorder=4)

ax_snp.legend(loc='best', frameon=True, fontsize=10, numpoints=1, ncol=1, handlelength=1)

if sweep_type == "partial":
	modification_difference_threshold = 50

# Now fill in the graphics
ax_snp.fill_between([1e-01,1], [ymin,ymin],[ymax,ymax],color='0.8',zorder=1)
ax_snp.fill_between([1e0,modification_difference_threshold],[ymin,ymin],[ymax,ymax],color='#deebf7',zorder=1)
ax_snp.fill_between([replacement_difference_threshold,1e05],[ymin,ymin],[ymax,ymax],color='#fee0d2',zorder=1)

ax_snp.text( exp((log(1e05)+log(replacement_difference_threshold))/2), ymax*1.2, 'putative\nreplacement',fontsize=12,fontstyle='italic',ha='center',color='#fc9272',zorder=1)
ax_snp.text( exp((log(1)+log(modification_difference_threshold))/2), ymax*1.2, 'putative\nmodification',fontsize=12,fontstyle='italic',ha='center',color='#9ecae1',zorder=1)

fig_snp.savefig('%s/temporal_snp_changes_%s_pooled.pdf' % (config.analysis_directory, sweep_type),bbox_inches='tight')