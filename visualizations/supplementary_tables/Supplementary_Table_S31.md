# Supplementary Table S31. Matched citrate contrasts in plate and flask endpoints

Citrate-dependent changes occur in the abiotic plate signal as well as the plate net response and apparent ferric generation in flasks. Each contrast uses zero citrate at the same gallium concentration. Plate endpoints are in AU and flask endpoints in mM.

**Figure panels** 7a–c · [All Supplementary tables](README.md)

**Column definitions**

- `readout` — Identifies apparent flask Fe(III) generation from Fe(II) depletion (mM), net biotic-minus-abiotic change in plate A380 (AU), or abiotic change in plate A380 (AU).
- `citric_acid_mM` — Citric acid concentration in the treated condition, in mM.
- `reference_citric_acid_mM` — Citric acid concentration of the matched reference, 0 mM for every contrast.
- `Ga_mM` — Ga(III) concentration in mM, held equal in the treated condition and its zero-citric-acid reference.
- `candidate_batch_n` — Number of experimental batches containing a treated vessel for this readout, citric acid concentration, and gallium concentration.
- `paired_batch_n` — Number of batches with at least one eligible treated–reference endpoint pair.
- `candidate_unit_n` — Number of treated physical microplates or shake flasks considered for the endpoint comparison.
- `paired_unit_n` — Number of treated units with a finite matched reference at the same final acquisition and elapsed time. Plate conditions are paired within physical plate and flask cultures by batch and starter culture.
- `unavailable_unit_n` — Number of candidate treated units lacking an eligible pair because the reference, final acquisition, timing, or measured endpoint could not be matched.
- `mean_effect` — Mean treated endpoint minus matched zero-citric-acid endpoint, in the units of the readout. Pair differences are averaged within batches, then batch means receive equal weight.
- `between_batch_effect_sd` — Sample SD of the batch mean endpoint differences, in the units of the readout, for two or more eligible batches.
- `within_batch_pair_sd` — Sample SD among unit-pair endpoint differences, in the units of the readout, when a single batch contributes at least two eligible pairs.
- `ci_low` — Lower limit of the two-sided 95% Student’s t interval for the mean endpoint effect, calculated from two or more independent batch means.
- `ci_high` — Upper limit of the same two-sided 95% Student’s t interval, in AU for plate readouts or mM for flasks.
- `endpoint_min_d` — Earliest elapsed day of a matched final acquisition among eligible treated–reference pairs.
- `endpoint_max_d` — Latest elapsed day of a matched final acquisition among eligible treated–reference pairs.
- `status` — Indicates whether a batch-based confidence interval was estimable, only one eligible batch contributed a descriptive mean, or no eligible endpoint pair was available.

Blank numeric cells indicate quantities that were unavailable or inapplicable for that row.

| readout | citric&#95;<wbr>acid&#95;<wbr>mM | reference&#95;<wbr>citric&#95;<wbr>acid&#95;<wbr>mM | Ga&#95;<wbr>mM | candidate&#95;<wbr>batch&#95;<wbr>n | paired&#95;<wbr>batch&#95;<wbr>n | candidate&#95;<wbr>unit&#95;<wbr>n | paired&#95;<wbr>unit&#95;<wbr>n | unavailable&#95;<wbr>unit&#95;<wbr>n | mean&#95;<wbr>effect | between&#95;<wbr>batch&#95;<wbr>effect&#95;<wbr>sd | within&#95;<wbr>batch&#95;<wbr>pair&#95;<wbr>sd | ci&#95;<wbr>low | ci&#95;<wbr>high | endpoint&#95;<wbr>min&#95;<wbr>d | endpoint&#95;<wbr>max&#95;<wbr>d | status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| flask apparent Fe(III) (mM) | 1.0 | 0.0 | 0.0 | 1 | 1 | 3 | 3 | 0 | 0.6333333333333352 |  | 1.5502687938978013 |  |  | 7.4965 | 7.4965 | single&#95;<wbr>batch&#95;<wbr>descriptive |
| flask apparent Fe(III) (mM) | 6.0 | 0.0 | 0.0 | 1 | 1 | 3 | 3 | 0 | -0.3333333333333333 |  | 0.37859388972001295 |  |  | 7.4965 | 7.4965 | single&#95;<wbr>batch&#95;<wbr>descriptive |
| flask apparent Fe(III) (mM) | 12.0 | 0.0 | 0.0 | 1 | 1 | 3 | 3 | 0 | 0.0 |  | 0.0 |  |  | 8.6255 | 8.6255 | single&#95;<wbr>batch&#95;<wbr>descriptive |
| flask apparent Fe(III) (mM) | 18.0 | 0.0 | 0.0 | 1 | 1 | 3 | 3 | 0 | 0.06666666666666288 |  | 0.05773502691895929 |  |  | 8.6255 | 8.6255 | single&#95;<wbr>batch&#95;<wbr>descriptive |
| flask apparent Fe(III) (mM) | 1.0 | 0.0 | 1.0 | 1 | 1 | 3 | 3 | 0 | -1.0 |  | 1.5620499351813284 |  |  | 7.4965 | 7.4965 | single&#95;<wbr>batch&#95;<wbr>descriptive |
| flask apparent Fe(III) (mM) | 6.0 | 0.0 | 1.0 | 1 | 1 | 3 | 3 | 0 | 1.1333333333333353 |  | 1.3316656236958848 |  |  | 7.4965 | 7.4965 | single&#95;<wbr>batch&#95;<wbr>descriptive |
| flask apparent Fe(III) (mM) | 12.0 | 0.0 | 1.0 | 1 | 1 | 3 | 3 | 0 | 0.0 |  | 0.0 |  |  | 8.6255 | 8.6255 | single&#95;<wbr>batch&#95;<wbr>descriptive |
| flask apparent Fe(III) (mM) | 18.0 | 0.0 | 1.0 | 1 | 1 | 3 | 3 | 0 | 0.06666666666666288 |  | 0.05773502691895929 |  |  | 8.6255 | 8.6255 | single&#95;<wbr>batch&#95;<wbr>descriptive |
| plate net A380 (AU) | 1.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | 0.14589999999999984 | 0.567419892535012 |  | -1.263649153443543 | 1.5554491534435428 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 3.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | -0.5253333333333337 | 0.07121463644816609 |  | -0.7022402973628126 | -0.3484263693038547 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 6.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | -0.7712111111111115 | 0.7176652920975641 |  | -2.553990527624987 | 1.0115683054027638 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 12.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | -1.192211111111111 | 0.8366810425945109 |  | -3.2706420417267195 | 0.8862198195044975 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 18.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | -1.4170444444444443 | 0.8754605266393819 |  | -3.5918089538181377 | 0.7577200649292488 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 24.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | -1.9255444444444445 | 0.28721716586995716 |  | -2.639031437644055 | -1.212057451244834 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 30.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | -2.0934000000000004 | 0.2787905805040367 |  | -2.785954194710844 | -1.4008458052891564 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 36.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | -2.2427000000000006 | 0.32114508836003935 |  | -3.0404686249385664 | -1.4449313750614348 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 1.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | -0.2157444444444443 | 0.02479539234376957 |  | -0.27733961364324766 | -0.1541492752456409 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 3.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | -0.5314777777777778 | 0.05230897876329884 |  | -0.6614204845868357 | -0.40153507096872004 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 6.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | -0.9788444444444445 | 0.20214434978041831 |  | -1.4809988469512312 | -0.47669004193765785 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 12.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | -1.422722222222222 | 0.29793342718976157 |  | -2.1628298842953297 | -0.6826145601491141 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 18.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | -1.724622222222222 | 0.33034323283742495 |  | -2.5452403047351893 | -0.9040041397092548 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 24.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | -1.687788888888889 | 0.2817387133366524 |  | -2.387666651548483 | -0.9879111262292951 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 30.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | -1.9678888888888888 | 0.27360966181229185 |  | -2.6475729680960574 | -1.2882048096817202 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate net A380 (AU) | 36.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | -2.188288888888889 | 0.27078250873539994 |  | -2.860949930520859 | -1.5156278472569191 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 1.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | -0.18721111111111108 | 0.5233918877362813 |  | -1.4873886374610026 | 1.1129664152387806 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 3.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | 0.42933333333333334 | 0.05284288031513802 |  | 0.29806434154498984 | 0.5606023251216768 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 6.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | 0.48461111111111127 | 0.6525505194008893 |  | -1.1364142429549031 | 2.1056364651771258 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 12.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | 0.9375444444444444 | 0.7219312342855356 |  | -0.8558321599347175 | 2.7309210488236064 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 18.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | 1.228511111111111 | 0.7307954375072967 |  | -0.5868853947758463 | 3.0439076169980686 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 24.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | 1.6872444444444445 | 0.14241465410449744 |  | 1.3334668314775835 | 2.041022057411306 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 30.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | 1.8950111111111114 | 0.14975092406813145 |  | 1.523009193264006 | 2.267013028958217 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 36.0 | 0.0 | 0.0 | 3 | 3 | 3 | 3 | 0 | 2.037577777777778 | 0.14259038978577487 |  | 1.6833636131777556 | 2.3917919423778002 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 1.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | 0.14176666666666668 | 0.0205979233041694 |  | 0.09059858860303852 | 0.19293474473029484 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 3.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | 0.4000333333333334 | 0.036852544009878044 |  | 0.30848653898445666 | 0.49158012768221016 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 6.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | 0.7472000000000002 | 0.11665022265064252 |  | 0.4574247828294665 | 1.036975217170534 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 12.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | 1.1451555555555553 | 0.18582267034679822 |  | 0.683546452448924 | 1.6067646586621867 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 18.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | 1.4447888888888887 | 0.21398159405096628 |  | 0.9132291414864332 | 1.9763486362913443 | 6.841712962962963 | 6.878344907407407 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 24.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | 1.5613888888888894 | 0.14116253962221229 |  | 1.2107217007269018 | 1.912056077050877 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 30.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | 1.7896555555555553 | 0.19646206683149722 |  | 1.3016167264110194 | 2.2776943847000912 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
| plate abiotic A380 (AU) | 36.0 | 0.0 | 1.0 | 3 | 3 | 3 | 3 | 0 | 1.8969222222222222 | 0.11592689906111502 |  | 1.6089438404582324 | 2.184900603986212 | 6.866412037037037 | 6.998194444444445 | marginal&#95;<wbr>paired&#95;<wbr>batch&#95;<wbr>t&#95;<wbr>interval |
