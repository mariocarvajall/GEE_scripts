// var pkg_trend = require('users/kongdd/public:Math/pkg_trend.js');
var pkg_trend = {};

/**
 * Add season, year, month, yearmonth, year-ingrow, ingrow,
 * 
 * add seasonal variables into img before regression
 * @param {[type]} img [description]
 * @param {boolean} pheno If true, 4-10 as growing season, spring:4-5, summer: 6-8, autumn:9-10
 *                        If false, just as traditional seasons.
 */
pkg_trend.addSeasonProb = function(img, pheno){
    if (pheno === undefined) {pheno = false;}
    var date  = ee.Date(img.get('system:time_start'));
    var month = date.get('month');
    var year  = date.get('year');
    var season;
    // year.subtract(1).multiply(10).add(4)
    var ingrow = ee.Algorithms.If(month.gte(4).and(month.lte(10)), "true", "false");

    if (pheno){
        /** 4-10 as growing season */
        season = ee.Algorithms.If(month.lte(3), ee.String(year.subtract(1)).cat("_winter"), season);
        season = ee.Algorithms.If(month.gte(4).and(month.lte(5)), ee.String(year).cat("_spring"), season);
        season = ee.Algorithms.If(month.gte(6).and(month.lte(8)), ee.String(year).cat("_summer"), season);
        season = ee.Algorithms.If(month.gte(9).and(month.lte(10)),ee.String(year).cat("_autumn"), season);
        season = ee.Algorithms.If(month.gte(11), ee.String(year).cat("_winter"), season);
    } else {
        /**traditional seasons*/
        season = ee.Algorithms.If(month.lte(2), ee.String(year.subtract(1)).cat("_winter"), season);
        season = ee.Algorithms.If(month.gte(3).and(month.lte(5)), ee.String(year).cat("_spring"), season);
        season = ee.Algorithms.If(month.gte(6).and(month.lte(8)), ee.String(year).cat("_summer"), season);
        season = ee.Algorithms.If(month.gte(9).and(month.lte(11)),ee.String(year).cat("_autumn"), season);
        season = ee.Algorithms.If(month.gte(12), ee.String(year).cat("_winter"), season);
    }
    
    return img.set('season', season)
        .set('ingrow', ingrow)
        .set('year-ingrow', year.format().cat('-').cat(ingrow))
        .set('year', year.format())
        .set('month', month.format("%02d"))
        .set('yearmonth', date.format('YYYY-MM')); //seasons.get(month.subtract(1))
}

/** add dn prop to every img */
pkg_trend.add_dn_date = function(img, beginDate, IncludeYear, n){
    beginDate = beginDate || img.get('system:time_start');
    if (IncludeYear === undefined) { IncludeYear = true; }
    n = n || 8;

    beginDate = ee.Date(beginDate);
    var year  = beginDate.get('year');
    // var month = beginDate.get('month');

    var diff  = beginDate.difference(ee.Date.fromYMD(year, 1, 1), 'day').add(1);
    var dn    = diff.subtract(1).divide(n).floor().add(1).int();
    
    var yearstr  = year.format('%d'); //ee.String(year);
    dn   = dn.format('%02d'); //ee.String(dn);
    dn   = ee.Algorithms.If(IncludeYear, yearstr.cat("-").cat(dn), dn);
    
    return ee.Image(img)
        .set('system:time_start', beginDate.millis())
        // .set('system:time_end', beginDate.advance(1, 'day').millis())
        .set('date', beginDate.format('yyyy-MM-dd')) // system:id
        .set('year', yearstr)
        .set('month', beginDate.format('MM'))
        .set('yearmonth', beginDate.format('YYYY-MM'))
        .set('dn', dn); //add dn for aggregated into 8days
}

/**
 * return a function used to add dn property
 * @param {boolean} IncludeYear [description]
 */
pkg_trend.add_dn = function(IncludeYear, n) {
    if (typeof IncludeYear === 'undefined') { IncludeYear = true; }
    if (typeof n === 'undefined') { n = 8; }
    return function (img) {
        return pkg_trend.add_dn_date(img, img.get('system:time_start'), IncludeYear, n);
    };
}

pkg_trend.imgcol_addSeasonProb = function (imgcol) {
    return imgcol.map(function (img) { return addSeasonProb(img, false); });
}

pkg_trend.imgcol_last = function(imgcol, n) {
    n = n || 1;
    // ee.Image(imgcol_grace.reduce(ee.Reducer.last())); properties are missing
    var res = imgcol.toList(n, imgcol.size().subtract(n));
    if (n <= 1) { res = ee.Image(res.get(0)); }
    return res;
}

pkg_trend.showdata = function(ImgCol) {
    ImgCol.filter(filter_date).aside(print);
}

pkg_trend.linearTrend = function (ImgCol, robust) {
    // img should only have dependant band
    function createConstantBand(img) {
        // img = ee.Image(img);
        var year = ee.Image(ee.Number.parse(img.get('Year'))).toInt16();
        return img.addBands(ee.Image([1, year]));
    }

    if (robust === undefined) { robust = false; }
    ImgCol = ImgCol.map(createConstantBand)
        .select([1, 2, 0], ['constant', 'Year', 'y']);

    var FUN;
    if (robust) {
        FUN = ee.Reducer.robustLinearRegression(2, 1);
    } else {
        FUN = ee.Reducer.linearRegression(2, 1);
    }
    var n = ee.Number(ImgCol.size());
    var bandnames = ['offset', 'slope'];

    var regression = ImgCol.reduce(FUN);
    var coef = regression.select('coefficients').arrayProject([0]).arrayFlatten([bandnames]);
    // root mean square of the residuals of each dependent variable
    // actually, it is RMSE, not residuals
    var RMSE = regression.select('residuals').arrayFlatten([['RMSE']]);
    var offset = coef.select('offset');
    var slope = coef.select('slope');

    /** try to get the tval to juage regression significant level */
    var Sx = n.multiply(n.add(1)).multiply(n.subtract(1)).divide(12);

    var tval, formula = false;
    if (formula) {
        // solution1: statistical formula
        ImgCol = ImgCol.map(function (img) {
            var pred = img.select(['Year']).multiply(slope).add(offset).rename('pred');
            var re = img.expression('pow(y - pred, 2)', { NDVI: img.select('y'), pred: pred }).rename('re');
            return img.addBands(pred).addBands(re);
        });
        var Sy = ImgCol.select('re').sum();
        tval = slope.expression('slope/sqrt(Sy/(Sx*(n-2)))', { slope: slope, Sx: Sx, Sy: Sy, n: n }).rename('tval');
    } else {
        // solution2: lazy method
        var adj = n.divide(n.subtract(2)).sqrt();
        tval = RMSE.expression('slope/(b()*adj)*sqrt(Sx)', { slope: slope, Sx: Sx, adj: adj }).rename('tval');
    }
    return coef.addBands(tval);
}

pkg_trend.imgcol_trend = function (imgcol, band, robust) {
    if (band === undefined) { band = 0; }
    if (robust === undefined) { robust = true; }

    // add seasonal variable
    imgcol = imgcol.select(band).map(function (img) { return addSeasonProb(img); });
    var trend = pkg_trend.linearTrend(imgcol, robust); //ee.Image
    return trend;
}

pkg_trend.seq_len = function (n) {
    return Array(n).join().split(',').map(function (e, i) { return i; });
};

pkg_trend.seq = function (from, to, by) {
    by = by || 1;
    var res = [];
    for (var i = from; i <= to; i += by) { res.push(i); }
    return res;
};

pkg_trend.seq_date = function (date_begin, date_end, by) {
    var day_secs = 86400000
    by = by * day_secs || day_secs;
    return ee.List.sequence(date_begin.millis(), date_end.millis(), by)
        .map(function (x) { return ee.Date(x); });
}

pkg_trend.date_format = function(dates) {
    return dates.map(function(date) { return ee.Date(date).format('yyyy-MM-dd'); })
}

// by: in day
pkg_trend.seq_yeardate = function (year, by) {
    var date_begin = ee.Date.fromYMD(year, 1, 1);
    var date_end = ee.Date.fromYMD(year, 12, 31);

    var day_secs = 86400000; // milliseconds
    by = by * day_secs || day_secs;
    return ee.List.sequence(date_begin.millis(), date_end.millis(), by)
        .map(function (x) { return ee.Date(x); });
}

/**
 * get empty imgcol by `seq.Date()`
 * @param {*} date_begin 
 * @param {*} date_end 
 * @param {*} by 
 */
pkg_trend.seq_imgcol = function (date_begin, date_end, by) {
    var dates = pkg_trend.seq_date(date_begin, date_end, by);
    /** blank ImgCol used to select the nearest Imgs */
    return dates.map(function (date) {
        return pkg_trend.add_dn_date(ee.Image(0), date);
    });
}

pkg_trend.YearDn2date = function(x, n){
    x = ee.String(x);
    n = n || 8;
    // var year = ee.Number.parse(x.slice(0,4));
    var i   = ee.Number.parse(x.slice(5,7));
    var doy = i.subtract(1).multiply(n).add(1);
    var datestr = x.slice(0, 5).cat(doy);
    return ee.Date.parse('Y-D', datestr); 
};

pkg_trend.copyProperties = function(source, destination, properties) {
    // properties = properties || ['system:time_start']; // , 'system:id'
    properties = properties || destination.propertyNames();
    return source.copyProperties(destination)
        .copyProperties(destination, properties);
};

/**
 * bandList -> reducerList
 * 
 * - 1 -> N : e.g. Tair to Tavg, Tmax, Tmin
 * - N -> 1 : 
 * - N -> N : N reduce corresponding to N List
 * 
 * @param {*} reducerList 
 * @param {*} bandList 
 * @param {*} ImgCol 
 */
pkg_trend.check_aggregate = function(bandList, reducerList, ImgCol) {
    function check_list(x) {
        if (!Array.isArray(x)) x = [x];
        return x;
    }

    reducerList = reducerList || ['mean'];
    bandList   = bandList || ImgCol.first().bandNames();

    reducerList = check_list(reducerList);
    bandList   = check_list(bandList);

    if (bandList.length === 1 && reducerList.length > 1) {
        temp = bandList[0];
        bandList = [];
        reducerList.forEach(function (reducer, i) {
            bandList.push(temp);
        });
    }
    return { bandList: bandList, reducerList: reducerList };
}

/**
 * aggregate_prop
 *
 * @param  {[type]} ImgCol  [description]
 * @param  {[type]} prop    The value of "prop" should be string.
 * @param  {[type]} reducer [description]
 * @param  {boolean} delta  If true, reducer will be ignore, and return
 *                          Just deltaY = y_end - y_begin. (for dataset like GRACE)
 * @return {[type]}         [description]
 */
pkg_trend.aggregate_process = function (ImgCol, prop, prop_val, bandList, reducerList, delta) {
    var nreducer = reducerList.length;
    var imgcol = ImgCol.filterMetadata(prop, 'equals', prop_val).sort('system:time_start');

    var first = ee.Image(imgcol.first());
    var last  = pkg_trend.imgcol_last(imgcol);

    var ans = ee.Image([]);
    if (!delta) {
        for (var i = 0; i < nreducer; i++) {
            var bands = bandList[i];
            var reducer = reducerList[i];
            var img_new = imgcol.select(bands).reduce(reducer);
            ans = ans.addBands(img_new);
        }
    } else {
        ans = last.subtract(first);
    }
    return pkg_trend.copyProperties(ee.Image(ans), first);
}

pkg_trend.aggregate_process2 = function (imgcol, bandList, reducerList) {
    var nreducer = reducerList.lengh;
    var first = ee.Image(imgcol.first());

    var ans = ee.Image([]);
    for (var i = 0; i < nreducer; i++) {
        var bands = bandList[i];
        var reducer = reducerList[i];
        var img_new = imgcol.select(bands).reduce(reducer);
        ans = ans.addBands(img_new);
    }
    return ee.Image(pkg_trend.copyProperties(ee.Image(ans), first));
};

pkg_trend.aggregate_prop = function (ImgCol, prop, reducerList, bandList, delta) {
    if (delta === undefined) { delta = false; }
    var dates = ee.Dictionary(ImgCol.aggregate_histogram(prop)).keys();
    var options = pkg_trend.check_aggregate(bandList, reducerList, ImgCol);

    function process(prop_val) {
        return pkg_trend.aggregate_process(ImgCol, prop, prop_val, 
            options.bandList, options.reducerList, delta);
    }

    var ImgCol_new = dates.map(process);
    var bands = ee.List(options.bandList).flatten();

    var out = ee.ImageCollection(ImgCol_new);
    if (options.reducerList.length === 1) {
        out = out.select(ee.List.sequence(0, bands.length().subtract(1)), bands);
    }
    return out;
};

exports = pkg_trend;
