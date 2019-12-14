figure(1)
boxplot(log10(transpose([
  0.0097159098775147501, 0.0097159098816989586, 0.0097159098775153052, 0.0097159098775144725, 0.0097159098814706968, 0.009715909877514417, 0.009715909877514417, 0.0097159098775146946, 2.9610439955529699e-07, 0.009715909877514417;
  0.0097159098775144725, 0.0097159098775144725, 9.1202721147487509e-07, 2.9604612805655961e-07, 2.8240574267979213e-07, 0.0097159098775143615, 0.009715909877514417, 8.703469654980367e-07, 1.9184732324983855e-07, 0.0097159098775559949;
  7.8925656793460419e-07, 4.5167366014009502e-07, 8.953617184892515e-07, 4.8080288883589617e-07, 3.0264105166377675e-07, 9.2243452370910362e-07, 4.1576706727441959e-07, 7.8025902139344794e-07, 5.8010319892165541e-07, 0.009715909877514417;
  1.6780861555876214e-07, 9.9551411641973786e-07, 9.4385666987717443e-07, 9.9795140334757448e-07, 6.1869975481743111e-07, 2.1678732686769564e-07, 7.8961181243508705e-07, 5.093169421765964e-07, 6.7679145526744477e-07, 2.2223144152677676e-07;
  8.0763948961948273e-07, 3.3563194828944987e-07, 9.6696112988903238e-07, 8.9223717525399593e-07, 4.8668165059106983e-07, 6.282652376121689e-07, 5.1246668403281959e-07, 7.1222727693331933e-07, 6.5849398855899466e-07, 8.5210616213027279e-07
])))
set(gca(), 'xtick', [1, 2, 3, 4, 5], 'xticklabel', {'50', '100', '250', '500', '1000'})
xlabel('Population size')
ylabel('Optimum unit fitness (log_{10})')
title('Optimum unit fitness (log_{10}) for different population sizes (mutation chance 50%)')

figure(2)
boxplot(log10(transpose([
  4.1521116533882463e-07, 0.0097159098775143615, 5.3244342113067944e-07, 9.7541401272716044e-07, 4.2986784698495484e-07, 9.2564973253050908e-07, 3.4199133952528271e-07, 9.8672833231949753e-07, 7.7055424307914677e-07, 2.7003221708676861e-07;
  1.068967826478584e-07, 5.2290610619287747e-08, 5.6807202053477113e-07, 1.7947002067808526e-07, 3.8555463510725474e-07, 3.7808732772370846e-07, 3.9903453785683496e-07, 5.7885673554469719e-07, 1.6960022491963045e-07, 4.288489354742353e-07;
  4.6911127110638162e-07, 9.3320959321241403e-07, 9.9403508080353475e-07, 7.2929035610513893e-07, 6.3511439724051044e-07, 8.7862347064593393e-07, 3.5249546515014885e-07, 3.4578690205622209e-07, 9.5206124017810367e-07, 2.9264487988456267e-07;
  8.5715658665330707e-07, 6.7469234610317841e-07, 3.9108502714046978e-08, 5.1772799330995767e-07, 7.9526764240611669e-07, 2.8741089674877429e-08, 8.642032399697186e-07, 3.2769314567415719e-07, 2.7369687627398775e-07, 1.9939900269827504e-07;
  8.6075318517231381e-07, 8.5507175878740682e-07, 6.5396149040441998e-07, 5.0856666283793928e-07, 1.0042627351936062e-07, 7.245343250916747e-07, 3.888474844893075e-07, 1.6887776543850208e-07, 8.6700553264540403e-07, 9.4792722166170407e-07;
  9.8082877553196113e-07, 7.17359726232214e-07, 1.9980492438342878e-07, 8.5094245910743993e-07, 7.873901747834644e-07, 5.4950546735099692e-07, 3.5739393922096241e-07, 3.1510020415126405e-07, 6.2964776392782085e-08, 4.1636499159558582e-07;
  2.0469542694190324e-07, 9.8093109884844765e-07, 5.3961287210801956e-07, 2.4436215684264795e-07, 9.2234830945114865e-07, 1.9741162776698573e-07, 4.8534457403048847e-07, 0.0097159098775143615, 4.2300488356517008e-07, 1.0528077187821339e-08;
  9.7405575677500522e-07, 1.8757938730074031e-07, 3.1356627483436483e-07, 1.8593546718959075e-07, 9.628617944623663e-07, 3.4372149743111535e-07, 3.5987257746006307e-07, 5.8143644154906582e-07, 4.3203273791814212e-08, 5.1123409799957997e-07;
  7.612836860904082e-09, 2.6118062851354651e-07, 9.4762436364526081e-07, 8.7306815593768405e-07, 3.652328774417235e-07, 8.8292947336698901e-07, 4.7649103052105346e-07, 3.8299960403165301e-07, 5.9281256470367083e-07, 5.3069582434117635e-07
])))
set(gca(), 'xtick', [1, 2, 3, 4, 5, 6, 7, 8, 9], 'xticklabel', {'10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%'})
xlabel('Mutation chance')
ylabel('Optimum unit fitness (log_{10})')
title('Optimum unit fitness (log_{10}) for different mutation chances (population size 250)')