// find $(pwd)/results -type f \( -name "*.pdf" -o -name "*.png" \) | sed "s|$(pwd)/||" | sed 's/.*/"&",/' >> slide-show.typ 
// typst compile slide-show.typ 

#let slide(body) = {
  pagebreak(weak: true)
  align(center + horizon)[#body]
}

#let scale = 150
#set page(height: scale*9pt, width: scale*16pt)

#let files = (
"results/005_1/component-graphs/current/sub_pathway-8.pdf",
"results/005_1/component-graphs/current/sub_pathway-9.pdf",
"results/005_1/component-graphs/current/chemical_name-8.pdf",
"results/005_1/component-graphs/current/chemical_name-9.pdf",
"results/005_1/component-graphs/current/chemical_name-7.pdf",
"results/005_1/component-graphs/current/sub_pathway-2.pdf",
"results/005_1/component-graphs/current/sub_pathway-3.pdf",
"results/005_1/component-graphs/current/chemical_name-6.pdf",
"results/005_1/component-graphs/current/chemical_name-4.pdf",
"results/005_1/component-graphs/current/sub_pathway-1.pdf",
"results/005_1/component-graphs/current/sub_pathway-0.pdf",
"results/005_1/component-graphs/current/chemical_name-5.pdf",
"results/005_1/component-graphs/current/chemical_name-1.pdf",
"results/005_1/component-graphs/current/sub_pathway-4.pdf",
"results/005_1/component-graphs/current/sub_pathway-5.pdf",
"results/005_1/component-graphs/current/chemical_name-0.pdf",
"results/005_1/component-graphs/current/chemical_name-2.pdf",
"results/005_1/component-graphs/current/sub_pathway-7.pdf",
"results/005_1/component-graphs/current/sub_pathway-6.pdf",
"results/005_1/component-graphs/current/chemical_name-3.pdf",
"results/005_1/component-graphs/former/sub_pathway-8.pdf",
"results/005_1/component-graphs/former/sub_pathway-9.pdf",
"results/005_1/component-graphs/former/chemical_name-8.pdf",
"results/005_1/component-graphs/former/chemical_name-9.pdf",
"results/005_1/component-graphs/former/chemical_name-7.pdf",
"results/005_1/component-graphs/former/sub_pathway-2.pdf",
"results/005_1/component-graphs/former/sub_pathway-3.pdf",
"results/005_1/component-graphs/former/chemical_name-6.pdf",
"results/005_1/component-graphs/former/chemical_name-4.pdf",
"results/005_1/component-graphs/former/sub_pathway-1.pdf",
"results/005_1/component-graphs/former/sub_pathway-0.pdf",
"results/005_1/component-graphs/former/chemical_name-5.pdf",
"results/005_1/component-graphs/former/chemical_name-1.pdf",
"results/005_1/component-graphs/former/sub_pathway-4.pdf",
"results/005_1/component-graphs/former/sub_pathway-5.pdf",
"results/005_1/component-graphs/former/chemical_name-0.pdf",
"results/005_1/component-graphs/former/chemical_name-2.pdf",
"results/005_1/component-graphs/former/sub_pathway-7.pdf",
"results/005_1/component-graphs/former/sub_pathway-6.pdf",
"results/005_1/component-graphs/former/chemical_name-3.pdf",
"results/005_1/chemical_name-only-curr-graph.pdf",
"results/005_1/chemical_name-both-graph.pdf",
"results/005_1/sub_pathway-both-graph.pdf",
"results/005_1/sub_pathway-only-form-graph.pdf",
"results/005_1/chemical_name-only-form-graph.pdf",
"results/005_1/sub_pathway-only-curr-graph.pdf",
"results/007/comp-metab-cor-hms/Former-10.pdf",
"results/007/comp-metab-cor-hms/Former-8.pdf",
"results/007/comp-metab-cor-hms/Former-9.pdf",
"results/007/comp-metab-cor-hms/Current-2.pdf",
"results/007/comp-metab-cor-hms/Current-3.pdf",
"results/007/comp-metab-cor-hms/Current-1.pdf",
"results/007/comp-metab-cor-hms/Current-4.pdf",
"results/007/comp-metab-cor-hms/Current-10.pdf",
"results/007/comp-metab-cor-hms/Current-5.pdf",
"results/007/comp-metab-cor-hms/Current-7.pdf",
"results/007/comp-metab-cor-hms/Current-6.pdf",
"results/007/comp-metab-cor-hms/Current-8.pdf",
"results/007/comp-metab-cor-hms/Current-9.pdf",
"results/007/comp-metab-cor-hms/Former-4.pdf",
"results/007/comp-metab-cor-hms/Former-5.pdf",
"results/007/comp-metab-cor-hms/Former-7.pdf",
"results/007/comp-metab-cor-hms/Former-6.pdf",
"results/007/comp-metab-cor-hms/Former-2.pdf",
"results/007/comp-metab-cor-hms/Former-3.pdf",
"results/007/comp-metab-cor-hms/Former-1.pdf",
"results/003/crpd-AIC-grid.png",
"results/003/log-AIC-grid.png",
"results/003/sqrt-AIC-grid.png",
"results/003/crpd-log-AIC-grid.png",
"results/003/crpd-lambda-grid-scatter.png",
"results/003/lambda-grid-scatter.png",
"results/003/crpd-sqrt-AIC-grid.png",
"results/003/AIC-grid.png",
"results/004/crpd-former-current-weight-scatter.png",
"results/004/crpd-log-partial-corr-distr.pdf",
"results/004/crpd-shortest-path-densities.pdf",
"results/004/shortest-path-densities.pdf",
"results/004/crpd-form-curr-metab-upset.pdf",
"results/004/former-current-weight-scatter.png",
"results/004/log-former-current-weight-scatter.png",
"results/004/crpd-partial-corr-distr.pdf",
"results/004/log-partial-corr-distr.pdf",
"results/004/crpd-log-former-current-weight-scatter.png",
"results/004/crpd-degree-densities.pdf",
"results/004/partial-corr-distr.pdf",
"results/004/degree-densities.pdf",
"results/005/form-curr-metab-upset.pdf"
)


#for ( f ) in files{
    slide()[
        #image(f, height:130%, fit:"contain")
    ]
}
