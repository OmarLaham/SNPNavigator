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

  //init genes and pathways datatable
  $('#dtGenesAndPathways').DataTable();


  function initSNPQuery() {

    //start loading spinner. Will need this only for the first query.
    $('#queryResultsSpinner').removeClass("d-none");

    let run_id = "run_1";
    let spec_chr = $("#lstSpecChr").val()
    var open_peak_cell_types = Array();
    let cpg_island = ($("#snpCpGIsland").is(':checked')) ? 1: 0;
    let close_to_another_ocr = ($("#snpCloseToAnotherOCR").is(':checked')) ? 1 : 0;
    var diseases_peaks_match = Array();
    var diseases_peaks_mismatch = Array();

    $('.atacCellType:checked').each(function() {
        open_peak_cell_types.push($(this).data('atac-cell-type'));
    });
    if(open_peak_cell_types.length > 0) {
      open_peak_cell_types = open_peak_cell_types.join(",")
    } else {
      open_peak_cell_types = "NA"
    }

    $('.diseaseMatch:checked').each(function() {
        diseases_peaks_match.push($(this).data('disease-match'));
    });
    if(diseases_peaks_match.length > 0) {
      diseases_peaks_match = diseases_peaks_match.join(",")
    } else {
      diseases_peaks_match = "NA"
    }

    $('.diseaseMismatch:checked').each(function() {
        diseases_peaks_mismatch.push($(this).data('disease-mismatch'));
    });
    if(diseases_peaks_mismatch.length > 0) {
      diseases_peaks_mismatch = diseases_peaks_mismatch.join(",")
    } else {
      diseases_peaks_mismatch = "NA"
    }

    let get_request_link = `/json_snp_query/${run_id}/${spec_chr}/${open_peak_cell_types}/${cpg_island}/${close_to_another_ocr}/${diseases_peaks_match}/${diseases_peaks_mismatch}`;

    $.get(get_request_link)
      .done(function(data, textStatus, jqXHR) {

        console.log(data); //TODO: remove this line =)

        let snps = data["snps"];

        //for consecutive queries, we must destroy the table before re-init
        if(tblSNPs) {
          tblSNPs.destroy();
        }
        tblSNPs = $('#dtSNPs').DataTable({data: snps});

        $("#queryResultsInstructions").hide(function() {
          $("#queryResultsWrapper").removeClass("d-none");
          $("#queryResultsWrapper").fadeIn(function() {});
        });


      })
      .fail(function(jqXHR, textStatus, errorThrown) {
        alert("Unable to get query results. Please try again or contact help.");
      });

  }

  //send initial SNP query
  $(document).ready(function() {

    $('#btnFilter').click(function() {
      initSNPQuery();
    });

  });
