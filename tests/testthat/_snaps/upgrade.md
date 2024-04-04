# Object from version 0.4 can be upgraded

    Code
      sim <- upgradeSim(simc.0.4)
    Condition
      Warning in `validSpeciesParams()`:
      The species parameter data frame is missing a `w_max` column. I am copying over the values from the `w_inf` column. But note that `w_max` should be the maximum size of the largest individual, not the asymptotic size of an average indivdidual.
    Message
      Initial effort has been set to 0.
    Condition
      Warning in `validSpeciesParams()`:
      The species parameter data frame is missing a `w_max` column. I am copying over the values from the `w_inf` column. But note that `w_max` should be the maximum size of the largest individual, not the asymptotic size of an average indivdidual.
      Warning in `validSpeciesParams()`:
      The species parameter data frame is missing a `w_max` column. I am copying over the values from the `w_inf` column. But note that `w_max` should be the maximum size of the largest individual, not the asymptotic size of an average indivdidual.

# Object from version 1.0 can be upgraded

    Code
      sim <- upgradeSim(simc.1.0)
    Condition
      Warning in `validSpeciesParams()`:
      The species parameter data frame is missing a `w_max` column. I am copying over the values from the `w_inf` column. But note that `w_max` should be the maximum size of the largest individual, not the asymptotic size of an average indivdidual.
    Message
      Initial effort has been set to 0.
    Condition
      Warning in `validSpeciesParams()`:
      The species parameter data frame is missing a `w_max` column. I am copying over the values from the `w_inf` column. But note that `w_max` should be the maximum size of the largest individual, not the asymptotic size of an average indivdidual.
      Warning in `validSpeciesParams()`:
      The species parameter data frame is missing a `w_max` column. I am copying over the values from the `w_inf` column. But note that `w_max` should be the maximum size of the largest individual, not the asymptotic size of an average indivdidual.

# Some functions work with params from earlier versions

    Code
      getEGrowth(params.0.4)
    Condition
      Warning in `validSpeciesParams()`:
      The species parameter data frame is missing a `w_max` column. I am copying over the values from the `w_inf` column. But note that `w_max` should be the maximum size of the largest individual, not the asymptotic size of an average indivdidual.
    Message
      Initial effort has been set to 0.
    Condition
      Warning in `validSpeciesParams()`:
      The species parameter data frame is missing a `w_max` column. I am copying over the values from the `w_inf` column. But note that `w_max` should be the maximum size of the largest individual, not the asymptotic size of an average indivdidual.
      Warning in `validSpeciesParams()`:
      The species parameter data frame is missing a `w_max` column. I am copying over the values from the `w_inf` column. But note that `w_max` should be the maximum size of the largest individual, not the asymptotic size of an average indivdidual.
      Warning in `validParams()`:
      Your MizerParams object was created with an earlier version of mizer. You can upgrade it with `params <- upgradeParams(params)` where you should replace `params` by the name of the variable that holds your MizerParams object.
    Output
                 w
      sp                1e-04    0.000129    0.000167    0.000215    0.000278
        Community 0.001901872 0.002274883 0.002721052 0.003254727 0.003893071
                 w
      sp             0.000359    0.000464    0.000599    0.000774       0.001
        Community 0.004656613 0.005569906 0.006662323 0.007968993 0.009531938
                 w
      sp             0.00129    0.00167    0.00215    0.00278    0.00359    0.00464
        Community 0.01140142 0.01363756 0.01631228 0.01951158 0.02333835 0.02791566
                 w
      sp             0.00599    0.00774       0.01     0.0129     0.0167     0.0215
        Community 0.03339071 0.03993958 0.04777286 0.05714247 0.06834973 0.08175505
                 w
      sp              0.0278    0.0359    0.0464  0.0599    0.0774       0.1
        Community 0.09778953 0.1169688 0.1399097 0.16735 0.2001721 0.2394315
                 w
      sp              0.129     0.167     0.215     0.278     0.359     0.464
        Community 0.2863908 0.3425601 0.4097459 0.4901086 0.5862328 0.7012097
                 w
      sp              0.599    0.774   1     1.29     1.67     2.15     2.78     3.59
        Community 0.8387368 1.003237 1.2 1.435354 1.716868 2.053594 2.456362 2.938124
                 w
      sp              4.64     5.99     7.74       10     12.9     16.7     21.5
        Community 3.514373 4.203642 5.028095 6.014247 7.193811 8.604721 10.29235
                 w
      sp              27.8    35.9     46.4     59.9     77.4      100      129
        Community 12.31097 14.7255 17.61359 21.06812 25.20017 30.14264 36.05446
                 w
      sp               167      215      278      359      464      599   774    1000
        Community 43.12576 51.58395 61.70102 73.80234 88.27707 105.5907 126.3 151.071
                 w
      sp              1290     1670     2150     2780     3590     4640     5990
        Community 180.7004 216.1408 258.5322 309.2376 369.8879 442.4334 529.2071
                 w
      sp              7740    10000   12900   16700   21500    27800   35900    46400
        Community 632.9996 757.1487 905.647 1083.27 1295.73 1549.859 1853.83 2217.418
                 w
      sp             59900    77400    1e+05   129000   167000   215000   278000
        Community 2652.316 3172.509 3794.726 4538.975 5429.189 6493.992 7767.619
                 w
      sp            359000   464000   599000   774000    1e+06  1290000  1670000
        Community 9291.013 11113.13 13292.51 15899.08 19016.37 22743.97 27200.27
                 w
      sp          2150000  2780000  3590000  4640000 5990000  7740000    1e+07
        Community 32525.3 38882.78 46459.78 55459.15 66073.5 78413.01 92324.49

