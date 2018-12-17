demandtime = ([0, 3600, 7200, 10800, 14400, 18000, 21600, 25200, 28800, 32400] - 3600)./3600;

demanddata = [1, 1, .9, .8, .7, .6, .5, .75, 1, 1];

demand = timeseries(demanddata,demandtime);