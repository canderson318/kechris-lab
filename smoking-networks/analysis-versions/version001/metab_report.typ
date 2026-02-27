// find $(pwd)/results/006/enrichment/ -type f | sed 's/.*/"&",/'  >> metab_report.typ 
// cd /Users/canderson/Documents/school/local-kechris-lab/kechris-lab/smoking-networks/analysis-versions/version001 && typst compile metab_report.typ
// cd /Users/canderson/Documents/school/local-kechris-lab/kechris-lab/smoking-networks/analysis-versions/version001 && open main.pdf


#set page(height: 11in, width: 8.5in)
#set text(size: 9pt)

#let is-int-like(s) = s.clusters().all(c => "0123456789".contains(c))
#let paths = read("metab_report_paths.txt").split("\n")       
#let today = datetime.today()

#set document(
  title: text(40pt)[Metabolite Reorganization Results],
  author: "Christian Anderson",
  date: today
)

#align(center + horizon)[
  #text(40pt, weight: "bold")[Metabolite Reorganization Results] \
  #v(2pt)
  Christian Anderson \
  #today.display("[day] [month repr:long] [year]")
]
#pagebreak()
#outline(indent: 1em)
#pagebreak(weak:true)

#heading(level: 1, strong[Node/Edge Overlaps])
#pagebreak(weak:true)

#let upset_pth = "results/005/node-upset.pdf"
#figure(image(upset_pth), caption:[Node intersection and difference])
#v(-1em)
#align(right, text(size: .5em, strong(upset_pth)))

#let upset_pth = "results/005/edge-upset.pdf"
#figure(image(upset_pth), caption:[Edge intersection and difference.])
#v(-1em)
#align(right, text(size: .5em, strong(upset_pth)))


#for ( path ) in paths {
    if path.ends-with(".csv"){
        let tab = csv(path)
        let num = path.split("/").last().split("-").first()
        
        if path.ends-with("superpathway.csv"){ 
            
            if is-int-like(num){
                pagebreak(weak:true)
                heading(level: 3)[Component #num]
            }else if num == "curr" {
                pagebreak(weak:true)
                heading(level: 3)[Unique To Current Smokers]
            }else if num == "form"{
                pagebreak(weak:true)
                heading(level: 3)[Unique To Former Smokers]
            }

            table(
                columns: 7,
                stroke: .5pt,
                inset: (x: 4pt, y: 5pt),
                table.hline(stroke: 1.5pt),
                [*Super-Pathway*], [*Expected Frequency*], [*Observed Frequency*],[*OR*],[*ln(OR)*], [*Fisher P*], [*FDR P*],
                table.hline(stroke: 0.75pt),
                ..tab.slice(1).flatten(),
                table.hline(stroke: 1.5pt),
            )

        } else{
            table(
                columns: 7,
                stroke: .5pt,
                inset: (x: 4pt, y: 5pt),
                table.hline(stroke: 1.5pt),
                [*Sub-Pathway*], [*Expected Frequency*], [*Observed Frequency*], [*OR*],[*ln(OR)*], [*Fisher P*], [*FDR P*],
                table.hline(stroke: 0.75pt),
                ..tab.slice(1).flatten(),
                table.hline(stroke: 1.5pt),
            )
        }
        v(-1em)
        align(right, text(size: .5em, strong[#path]))
        v(1em)
    } else if path.ends-with("pdf"){
        stack(
            align(center)[#image(path, width: 100%)],
            v(-5em),
            align(right, text(size: .5em, strong[#path])),
        )
        v(5em)
    } else if path.ends-with("txt"){
        let metabs = csv(path, delimiter: "\t")
        let chems = metabs.map(row => row.map(cell => cell.split("•").slice(-2).join("\n       ")))
        text[=== Unique Edges (in orange)]
        stack(
            table(columns: chems.at(0).len(), ..chems.flatten()),
            v(1em),
            align(right, text(size: .5em, strong[#path])),
        )
    }else{
        pagebreak(weak:true)
        heading(level:1)[#strong[#path]]
        v(1em)
    }
}


#pagebreak(weak:true)

#heading(level: 2,strong[Smoker Component Pairs])

#image("results/005_2/component-pairs-graph.pdf")

#heading(level: 2,strong[Current Smoker Component Metabolites])
#v(.5em)
#let lines = read("results/005/curr-metabs/comp-metabs.txt").split("\n")
#let i = 0
#for line in lines{
    i+=1
    let metabolites = line.split("||").map(chunk => chunk.split("•").last())
    if metabolites.len()!=1{
        text(strong[Component #i: ])
        for metab in metabolites{
            [#metab, ] 
        }
        v(.2em)
    }
}

#heading(level: 2,strong[Former Smoker Component Metabolites])
#v(.5em)
#let lines = read("results/005/form-metabs/comp-metabs.txt").split("\n")
#let i = 0
#for line in lines{
    i+=1
    let metabolites = line.split("||").map(chunk => chunk.split("•").last())
    if metabolites.len()!=1{
        text(strong[Component #i: ])
        for metab in metabolites{
            [#metab, ] 
        }
        v(.2em)
    }
}