Highcharts.chart('container', {

    title: {
        text: 'SNPs Dependency Wheel'
    },

    accessibility: {
        point: {
            valueDescriptionFormat: '{index}. From {point.from}-{point.fromsnp} to {point.to}-{point.tosnp}: {point.weight}.'
        }
    },

    series: [{
        keys: ['from', 'to', 'fromsnp', 'tosnp', 'weight'],
        data: [
        		['chr2', 'chr7', 'rs17480230', 'rs6949503', 1],
	       		['chr6', 'chr8', 'rs16901784', 'rs788001', 1],
	       		['chr7', 'chr11', 'rs66638610', 'rs6485671', 1],
            ['chr8', 'chr2', 'rs2056459', 'rs788001', 1],
            ['chr8', 'chr6', 'rs2056459', 'rs16901784', 2]
            ['chr10', 'chr2', 'rs17115213', 'rs788001',  1],
            ['chr13', 'chr13', 'rs2783125', 'rs2783126', 1]
        ],
        type: 'dependencywheel',
        name: 'SNPs Dependency Wheel',
        dataLabels: {
            color: '#333',
            style: {
                textOutline: 'none'
            },
            textPath: {
                enabled: true
            },
            distance: 10
        },
        tooltip: {
            pointFormatter: function() {
            	return `From <b>${this.from}</b>-${this.fromsnp} to <b>${this.to}</b>-${this.tosnp}: <b>${this.weight}</b>.`
            },
        },
        size: '95%'
    }]

});


document.getElementById('container').setAttribute("style","height:600px");
