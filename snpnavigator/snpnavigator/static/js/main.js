  /**
  * This function gets genomic coordinate data for ACMG genes from NCBI EUtils.
  * EUtils docs: https://www.ncbi.nlm.nih.gov/books/NBK25500/
  * Ideogram.js then draws those genes as annotations on chromosomes.
  */
  function getAndDrawAcmgGenes() {

    var acmgGenes, i, annots, geneClause, geneID, geneIDs, topWeight,
        eutils, esearch, esummary, url, defaultParams,
        ideo = this;

    ideo.annotDescriptions = {annots: {}}

    geneIDs = [];
    annots = [];

    // EUtils is a web API to access NCBI data
    eutils = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    defaultParams = '?db=gene&retmode=json&retmax=100';

    // Use ESearch to get NCBI Gene ID for gene name, e.g. BRCA1 -> 672
    esearch = eutils + 'esearch.fcgi' + defaultParams;

    // Use ESummary to get genomic coordinates for given gene IDs
    esummary = eutils + 'esummary.fcgi' + defaultParams;

    acmgGenes = [
      'APC', 'MYH11', 'ACTA2', 'TMEM43', 'DSP', 'PKP2', 'DSG2', 'DSC2',
      'BRCA1', 'BRCA2', 'SCN5A', 'RYR2', 'LMNA', 'MYBPC3', 'COL3A1', 'GLA',
      'APOB', 'LDLR', 'MYH7', 'TPM1', 'PRKAG2', 'TNNI3', 'MYL3', 'MYL2',
      'ACTC1', 'RET', 'PCSK9', 'BMPR1A', 'SMAD4', 'TNNT2', 'TP53', 'TGFBR1',
      'TGFBR2', 'SMAD3', 'KCNQ1', 'KCNH2', 'MLH1', 'MSH2', 'MSH6', 'PMS2',
      'RYR1', 'CACNA1S', 'FBN1', 'MEN1', 'MUTYH', 'NF2', 'OTC', 'SDHD',
      'SDHAF2', 'SDHC', 'SDHB', 'STK11', 'PTEN', 'RB1', 'TSC1', 'TSC2',
      'VHL', 'WT1', 'ATP7B'
    ];

    // Batch request all ACMG genes
    geneClause = '(' + acmgGenes.join('[symbol] OR ') + '[symbol])';

    topWeight = 0;
    weights = [];

    url = esearch + '&term=Homo-sapiens[Organism] AND ' + geneClause;
    d3.json(url).then(function(data) {
      geneIDs = data.esearchresult.idlist;

      // Batch request genomic coordinates for all ACMG genes
      url = esummary + '&id=' + geneIDs.join(',');
      d3.json(url).then(function(data) {
        var result, gene, chr, loc, start, stop, weight;
        for (i = 0; i < geneIDs.length; i++) {
          geneID = geneIDs[i];
          result = data.result[geneID];
          if (result.currentid !== '' || acmgGenes.indexOf(result.name) === -1 ) {
            // currentid case occurs when one gene symbol in a taxon maps to
            // multiple gene.  It seems to be annotation-run noise.
            // result.name case occurs when gene has a non-canonical alias
            // matching the gene symbol.
            // Ignore both.
            continue;
          }
          gene = result.name;
          chr = result.chromosome;
          loc = result.locationhist[0]; // better than 'genomicinfo'
          start = loc.chrstart;
          stop = loc.chrstop;
          weight = result.geneweight;

          if (result.weight > topWeight) {
            topWeight = result;
          }

          ////
          // This is the annotations API
          ////
          annots.push({
            name: gene, // required
            chr: chr, // required
            start: start, // required
            stop: stop, // required
            weight: weight // optional
          });

          var clinvarText = 'Pathogenic variants';
          var clinvarBase = 'https://www.ncbi.nlm.nih.gov/clinvar';

          var clinvarUrl =
            'https://www.ncbi.nlm.nih.gov/clinvar' +
            '?term=' + gene + '%5Bgene%5D%20AND%20%22clinsig%20pathogenic%22%5BProperties%5D';
          var clinvarLink =
            'ClinVar: <a target="_blank" href="' + clinvarUrl + '">' + clinvarText + '</a>';
          var fullName = result.description

          // Add properties like fullName or clinvarLink for later use in hyperlinkGene
          ideo.annotDescriptions.annots[gene] = {
            fullName,
            clinvarLink
          }

        }

        // NCBI assigns each gene a weight based on how well it is characterized.
        // http://www.ncbi.nlm.nih.gov/books/NBK3841/#EntrezGene.Sort_by
        // Gene weight is a rough proxy for general biomedical relevance.
        // Higher weight, more important.
        annots.sort(function(a, b) {
          return b.weight - a.weight;
        });

        // This block highlights the 10 most important ACMG genes.
        for (i = 0; i < annots.length; i++) {
          color = (i < 10) ? '#825' : '#BAB';
          annots[i].shape = 'triangle';
          annots[i].color = color;

          // BRCA1 is the most widely-recognized gene name.  Make it stand out.
          if (annots[i].name === 'BRCA1') {
            annots[i].shape = 'triangle';
            annots[i].color = '#F00';
          }
        }

        ideogram.drawAnnots(annots);

      });
    });
  }

  function enhanceTooltipContent(annot) {
    var term = '(' + annot.name + '[gene])+AND+(Homo+sapiens[orgn])';
    var url = 'https://ncbi.nlm.nih.gov/gene/?term=' + term;
    var description = this.annotDescriptions.annots[annot.name];

    annot.displayName =
      '<a target="_blank" href="' + url + '">' + annot.name + '</a><br/>' +
      description.fullName + '<br/><br/>' +
      description.clinvarLink + '<br/>';

    return annot
  }


  //init ideogram
  var d3 = Ideogram.d3;

  //datatable variables
  var tblSNPs;

  var config = {
    container: '#ideogramContainer',
    organism: 'human',
    resolution: 550,
    chrHeight: 350,
    annotationHeight: 4,
    onLoad: getAndDrawAcmgGenes,
    onWillShowAnnotTooltip: enhanceTooltipContent
  };
  var ideogram = new Ideogram(config);
  var pltManhattan = null;

  //init genes and pathways datatable
  $('#dtGenesAndPathways').DataTable();


  function plotManhattan(dataManhattan) {

    let peaks = dataManhattan["peaks"];
    let snpsDetails = dataManhattan["snps_details"];
    let series = dataManhattan["series"]; //selected and unselected SNPs series

    // Highcharts.setOptions({
    //     colors: ['rgba(5,141,199,0.5)', 'rgba(80,180,50,0.5)']
    // });


    // destroy if exist
    if(pltManhattan != null) {
        pltManhattan.destroy();
    }
    pltManhattan = Highcharts.chart('manhattanContainer', {
        chart: {
            type: 'scatter',
            zoomType: 'xy'
        },
        title: {
            text: '',
            align: 'left'
        },
        subtitle: {
            text:
            '',
            align: 'left'
        },
        xAxis: {
            visible: false,
            title: {
                text: 'DNA / Part of DNA'
            },
            labels: {
                format: ''
            },
            startOnTick: true,
            endOnTick: true,
            showLastLabel: true
        },
        yAxis: {
            title: {
                text: '-log10(p-val)'
            },
            labels: {
                format: '{value}'
            },
            min: 8,//using GWAS pval thresh
        },
        legend: {
            enabled: true
        },
        plotOptions: {
            scatter: {
                marker: {
                    radius: 2.5,
                    symbol: 'circle',
                    states: {
                        hover: {
                            enabled: true,
                            lineColor: 'rgb(100,100,100)'
                        }
                    }
                },
                states: {
                    hover: {
                        marker: {
                            enabled: false
                        }
                    }
                },
                jitter: {
                    x: 0.005
                }
            }
        },
        tooltip: {
            useHTML: true,
            style: {
              pointerEvents: 'auto'
            },
            formatter: function () {
                let snp_id = snpsDetails[this.x]["id"];
                let snp_chr = snpsDetails[this.x]["chr"];
                let snp_pos = snpsDetails[this.x]["pos"];
                return '<b>SNP ID</b>: ' + snp_id + '<br/>' +
                    '<b>Chr</b>: ' + snp_chr + '<br/>' +
                    '<b>Pos</b>: ' + snp_pos + '<br/>' +
                    '<b>-log10(p-val)</b>: ' + this.y + '<br />' +
                    '<hr />' +
                    '<b>dbSNP Search</b>: <a target="_blank" href="https://www.ncbi.nlm.nih.gov/snp/?term=' + snp_id + '">(Click)</a><br />' +
                    '<b>ClinVar</b>: <a target="_blank" href="https://www.ncbi.nlm.nih.gov/clinvar/?term=' + snp_id + '">(Click)</a><br />' +
                    '<b>GWAS Catalog</b>: <a target="_blank" href="https://www.ebi.ac.uk/gwas/variants/' + snp_id + '">(Click)</a><br />';

            }
        },
        series
    });


    return;
  }

  var open_peak_cell_types = "NA";
  let update_open_peak_cell_types = function() {
        open_peak_cell_types = Array();
        // update the status of overall peak cell_types
        $('.atacCellType').each(function() {
            if ($(this).prop('checked')) {
                open_peak_cell_types.push($(this).data('atac-cell-type'));
            }
        });
        if(open_peak_cell_types.length > 0) {
          open_peak_cell_types = open_peak_cell_types.join(",")
        } else {
          open_peak_cell_types = "NA"
        }
  }

  function initSNPQuery() {

    //start loading spinner and hide results from previous query (if available)
    $('#queryResultsSpinner').removeClass("d-none");
    $("#queryResultsWrapper").addClass("d-none");

    const current_url = window.location.href;
    const current_url_splitted = current_url.split("/");
    const run_id = current_url_splitted[current_url_splitted.length - 1];
    const pval_thresh = $('#lstPValThresh').val();
    const spec_chr = $("#lstSpecChr").val();
    const spec_gen_region = $('#lstGenomicRegion').val();
    const filter_eQTL = $('#lsteQTL').val();
    var cell_specific_ocrs = $('#lstCellSpecificOCRs').val();
    const cpg_island = ($("#snpCpGIsland").is(':checked')) ? 1: 0;
    const close_to_another_ocr = ($("#snpCloseToAnotherOCR").is(':checked')) ? 1 : 0;
    var condition_2_match = $("input:radio[name ='radioCondition2Match']:checked").val();
    var condition_3_match = "NA";// TODO: remove condition 3 if not needed //$("input:radio[name ='radioCondition3Match']:checked").val();
    let reverse_results = ($("#toggleReverseResults").is(':checked')) ? 1 : 0;

    let get_request_link = `/json_snp_query/${run_id}/${pval_thresh}/${spec_chr}/${spec_gen_region}/${filter_eQTL}/${open_peak_cell_types}/${cell_specific_ocrs}/${cpg_island}/${close_to_another_ocr}/${condition_2_match}/${condition_3_match}/${reverse_results}`;

    let TEST_MODE = false;
    if(TEST_MODE) {
        alert(get_request_link);
        return;
    }

    $.get(get_request_link)
      .done(function(data, textStatus, jqXHR) {

        console.log(data); //TODO: remove this line =)

        let selected_snps = data["selected_snps"];

        //for consecutive queries, we must destroy the table before re-init
        if(tblSNPs) {
          tblSNPs.destroy();
        }
        tblSNPs = $('#dtSNPs').DataTable({
            data: selected_snps,
            order: [[1, 'asc'], [2, 'asc']],
            dom: 'Bfrtip',
            buttons: ['csv', 'excel']
        });

        //plot Manhattan
        let dataManhattan = data["manhattan"];
        plotManhattan(dataManhattan);

        $("#queryResultsSpinner").addClass("d-none");
        $("#queryResultsWrapper").removeClass("d-none");
        $("#queryResultsWrapper").fadeIn(function() {});

      })
      .fail(function(jqXHR, textStatus, errorThrown) {
        alert("Unable to get query results. Please try again or contact help.");
      });

  }

  //send initial SNP query
  $(document).ready(function() {

    $('.atacCellType').change(function() {

        update_open_peak_cell_types();

        //disable other peak filters if no peak cell type is checked, and enable otherwise
        if (open_peak_cell_types == "NA") {
            $('#snpCpGIsland').attr("disabled", true);
            $('#snpCloseToAnotherOCR').attr("disabled", true);
        } else {
            $('#snpCpGIsland').attr("disabled", false);
            $('#snpCloseToAnotherOCR').attr("disabled", false);
        }

    });

    $('#btnFilter').click(function() {
      initSNPQuery();
    });

    //enable bootstrap toggles
    $(function () {
      $('[data-toggle="tooltip"]').tooltip();
    })


  });
