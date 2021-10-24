"""
solving systems of polynomial equations
"""

function test_pole28sys()


@polyvar x11 x12 x13 x14 x15 x16 x17 x18 x21 x22 x23 x24 x25 x26 x27 x28 # variables
x=[x11; x12; x13; x14; x15; x16; x17;x18; x21; x22; x23; x24; x25; x26; x27;x28]


#pole28sys

h=[4.02906133773525E-02*x11 + 1.03330704566256E-02*x12 - 5.24909824793066E-02*x13+
-4.20786318458456E-02*x14 + 4.57978985467796E-02*x15 + 8.42266523767165E-02*x16+
-2.93648884023350E-03*x17 - 1.07061117512882E-01*x18+
-4.45501027218261E-02*x21 + 5.26007926702697E-02*x22 + 1.36868129060802E-01*x23+
+2.88566286247173E-02*x24 - 1.81717394849268E-01*x25 - 1.84793108035429E-01*x26+
+9.84487400607885E-02*x27 + 3.13025292749779E-01*x28+
+7.88766907304129E-02*x11*x22+ 9.71113865576172E-02*x11*x23+
-2.17690682981990E-02*x11*x24-1.61480314075999E-01*x11*x25-1.12922435970517E-01*x11*x26+
 1.17283171143223E-01*x11*x27+ 2.39792426470830E-01*x11*x28-7.88766907304129E-02*x12*x21+
 1.27666802746855E-01*x12*x23+ 7.67941129360537E-02*x12*x24-1.31072071130368E-01*x12*x25+
-1.93850488813448E-01*x12*x26+ 3.58275953598466E-02*x12*x27+ 2.71090901265450E-01*x12*x28+
-9.71113865576172E-02*x13*x21-1.27666802746855E-01*x13*x22+ 1.29781942425538E-01*x13*x24+
 9.99925930608430E-02*x13*x25-5.58928544539397E-02*x13*x26-1.45719855980415E-01*x13*x27+
-5.43572386936603E-02*x13*x28+ 2.17690682981990E-02*x14*x21-7.67941129360537E-02*x14*x22+
-1.29781942425538E-01*x14*x23+ 1.93391155273753E-01*x14*x25+ 1.63441476970872E-01*x14*x26+
-1.24074557007834E-01*x14*x27-3.08279198824241E-01*x14*x28+ 1.61480314075999E-01*x15*x21+
 1.31072071130368E-01*x15*x22-9.99925930608430E-02*x15*x23-1.93391155273753E-01*x15*x24+
 2.09213395045194E-01*x15*x26+ 1.21545373042473E-01*x15*x27-1.56519673762409E-01*x15*x28+
 1.12922435970517E-01*x16*x21+ 1.93850488813448E-01*x16*x22+ 5.58928544539397E-02*x16*x23+
-1.63441476970872E-01*x16*x24-2.09213395045194E-01*x16*x25+ 2.36947830083076E-01*x16*x27+
 2.01220842276990E-01*x16*x28-1.17283171143223E-01*x17*x21-3.58275953598466E-02*x17*x22+
 1.45719855980415E-01*x17*x23+ 1.24074557007834E-01*x17*x24-1.21545373042473E-01*x17*x25+
-2.36947830083076E-01*x17*x26+ 2.94170740781183E-01*x17*x28-2.39792426470830E-01*x18*x21+
-2.71090901265450E-01*x18*x22+ 5.43572386936603E-02*x18*x23+ 3.08279198824241E-01*x18*x24+
 1.56519673762409E-01*x18*x25-2.01220842276990E-01*x18*x26-2.94170740781183E-01*x18*x27+
 + 3.27049414371538E-02;
 1.01563324725523E-01*x11*x22+ 8.45949419436556E-02*x11*x23+
-1.08531639783178E-01*x11*x24-1.90918944016144E-01*x11*x25+ 4.39987376158974E-02*x11*x26+
 2.77084122004202E-01*x11*x27+ 8.66475703014746E-02*x11*x28-1.01563324725523E-01*x12*x21+
 1.28854924474484E-01*x12*x23+ 2.65015340772530E-02*x12*x24-1.57243724827240E-01*x12*x25+
-8.24044359170837E-02*x12*x26+ 1.66656144182290E-01*x12*x27+ 1.61736579106925E-01*x12*x28+
-8.45949419436556E-02*x13*x21-1.28854924474484E-01*x13*x22+ 1.59769602146143E-01*x13*x24+
 1.11249039648114E-01*x13*x25-1.24458829197799E-01*x13*x26-2.12728234632579E-01*x13*x27+
 2.47838518424272E-02*x13*x28+ 1.08531639783178E-01*x14*x21-2.65015340772530E-02*x14*x22+
-1.59769602146143E-01*x14*x23+ 2.17849940042303E-01*x14*x25+ 7.65773918116980E-02*x14*x26
-2.50391752912112E-01*x14*x27-1.95442889784245E-01*x14*x28+ 1.90918944016144E-01*x15*x21+
 1.57243724827240E-01*x15*x22-1.11249039648114E-01*x15*x23-2.17849940042303E-01*x15*x24+
 2.23024338157211E-01*x15*x26+ 1.15710315757196E-01*x15*x27-1.69882093138171E-01*x15*x28+
-4.39987376158974E-02*x16*x21+ 8.24044359170837E-02*x16*x22+ 1.24458829197799E-01*x16*x23+
-7.65773918116980E-02*x16*x24-2.23024338157211E-01*x16*x25+ 2.97012930768203E-01*x16*x27+
 1.40369070231063E-01*x16*x28-2.77084122004202E-01*x17*x21-1.66656144182290E-01*x17*x22+
 2.12728234632579E-01*x17*x23+ 2.50391752912112E-01*x17*x24-1.15710315757196E-01*x17*x25+
-2.97012930768203E-01*x17*x26+ 2.99067484553945E-01*x17*x28-8.66475703014746E-02*x18*x21+
-1.61736579106925E-01*x18*x22-2.47838518424272E-02*x18*x23+ 1.95442889784245E-01*x18*x24+
 1.69882093138171E-01*x18*x25-1.40369070231063E-01*x18*x26-2.99067484553945E-01*x18*x27+
 4.05332206801165E-02*x11+ 2.97119117846005E-03*x12-4.89503209710238E-02*x13+
-1.37516252406762E-02*x14+ 5.71696321720076E-02*x15+ 3.41742046796462E-02*x16+
-5.84053386089786E-02*x17-6.20131132325031E-02*x18-3.19758471105613E-02*x21+
 9.05845150482069E-02*x22+ 1.16018623728405E-01*x23-8.84559163416066E-02*x24+
-2.19787027666882E-01*x25+ 1.32986259479277E-02*x26+ 2.99601281150956E-01*x27+
 1.28201713529575E-01*x28+ 3.70870932365451E-02;
 1.08709131095369E-01*x11*x22+ 3.91538123650870E-02*x11*x23+
-1.69364743140684E-01*x11*x24-1.01560045599980E-01*x11*x25+ 2.12972060199504E-01*x11*x26+
 1.83090657492785E-01*x11*x27-2.34034625818805E-01*x11*x28-1.08709131095369E-01*x12*x21+
 1.14026617666918E-01*x12*x23+ 2.14464645055560E-03*x12*x24-1.22306746548977E-01*x12*x25+
-7.35854331903625E-03*x12*x26+ 1.32709427652704E-01*x12*x27+ 1.65814269816207E-02*x12*x28+
-3.91538123650870E-02*x13*x21-1.14026617666918E-01*x13*x22+ 1.78421625687147E-01*x13*x24+
 6.24765649036638E-02*x13*x25-2.26039877781208E-01*x13*x26-1.44248493322658E-01*x13*x27+
 2.51454524605332E-01*x13*x28+ 1.69364743140684E-01*x14*x21-2.14464645055560E-03*x14*x22+
-1.78421625687147E-01*x14*x23+ 1.92552924432731E-01*x14*x25+ 7.26275721448962E-03*x14*x26+
-2.10368371313782E-01*x14*x27-2.12161349174562E-02*x14*x28+ 1.01560045599980E-01*x15*x21+
 1.22306746548977E-01*x15*x22-6.24765649036638E-02*x15*x23-1.92552924432731E-01*x15*x24+
 2.46485769078200E-01*x15*x26+ 8.20101037297266E-02*x15*x27-2.78799249324166E-01*x15*x28+
-2.12972060199504E-01*x16*x21+ 7.35854331903625E-03*x16*x22+ 2.26039877781208E-01*x16*x23+
-7.26275721448962E-03*x16*x24-2.46485769078200E-01*x16*x25+ 2.72384485564399E-01*x16*x27+
 1.66428221332583E-02*x16*x28-1.83090657492785E-01*x17*x21-1.32709427652704E-01*x17*x22+
 1.44248493322658E-01*x17*x23+ 2.10368371313782E-01*x17*x24-8.20101037297266E-02*x17*x25+
-2.72384485564399E-01*x17*x26+ 3.13630559530931E-01*x17*x28+ 2.34034625818805E-01*x18*x21+
-1.65814269816207E-02*x18*x22-2.51454524605332E-01*x18*x23+ 2.12161349174562E-02*x18*x24+
 2.78799249324166E-01*x18*x25-1.66428221332583E-02*x18*x26-3.13630559530931E-01*x18*x27+
 3.74377028051460E-02*x11+ 2.18745247212799E-04*x12-3.91901753863474E-02*x13+
-1.07937914564429E-03*x14+ 4.19161463730742E-02*x15+ 2.96270957710149E-03*x16+
-4.53346085200868E-02*x17-6.18130685742237E-03*x18-1.33405399362428E-02*x21+
 1.06877432788035E-01*x22+ 5.24871787561200E-02*x23-1.66247839790050E-01*x24+
-1.14858014764337E-01*x25+ 2.08480556143890E-01*x26+ 1.96291467375440E-01*x27+
-2.28056415798452E-01*x28+ 3.68337388465457E-02;
 1.05212222080466E-01*x11*x22-1.44885704440501E-02*x11*x23+
-1.73689723378655E-01*x11*x24+ 3.84833285387582E-02*x11*x25+ 2.39685079543829E-01*x11*x26+
-7.17561284334596E-02*x11*x27-3.02232987998254E-01*x11*x28-1.05212222080466E-01*x12*x21+
 1.05961012461766E-01*x12*x23-1.14698973593238E-04*x12*x24-1.07197141423625E-01*x12*x25+
 4.00302643023958E-04*x12*x26+ 1.08902924972498E-01*x12*x27-9.22893736977889E-04*x12*x28+
 1.44885704440501E-02*x13*x21-1.05961012461766E-01*x13*x22+ 1.74941659852867E-01*x13*x24+
-2.39953027322035E-02*x13*x25-2.41446031755296E-01*x13*x26+ 5.72700034270778E-02*x13*x27+
 3.04511055702801E-01*x13*x28+ 1.73689723378655E-01*x14*x21+ 1.14698973593238E-04*x14*x22+
-1.74941659852867E-01*x14*x23+ 1.77008483148401E-01*x14*x25-3.99543150919627E-04*x14*x26+
-1.79860750905986E-01*x14*x27+ 1.19407557306491E-03*x14*x28-3.84833285387582E-02*x15*x21+
 1.07197141423625E-01*x15*x22+ 2.39953027322035E-02*x15*x23-1.77008483148401E-01*x15*x24+
 2.44353363504262E-01*x15*x26-3.32765978886692E-02*x15*x27-3.08272439635527E-01*x15*x28+
-2.39685079543829E-01*x16*x21-4.00302643023958E-04*x16*x22+ 2.41446031755296E-01*x16*x23+
 3.99543150919627E-04*x16*x24-2.44353363504262E-01*x16*x25+ 2.48365920667172E-01*x16*x27+
-9.52543277497450E-04*x16*x28+ 7.17561284334596E-02*x17*x21-1.08902924972498E-01*x17*x22+
-5.72700034270778E-02*x17*x23+ 1.79860750905986E-01*x17*x24+ 3.32765978886692E-02*x17*x25+
-2.48365920667172E-01*x17*x26+ 3.13464339461298E-01*x17*x28+ 3.02232987998254E-01*x18*x21+
 9.22893736977889E-04*x18*x22-3.04511055702801E-01*x18*x23-1.19407557306491E-03*x18*x24+
 3.08272439635527E-01*x18*x25+ 9.52543277497450E-04*x18*x26-3.13464339461298E-01*x18*x27+
 3.52375344538437E-02*x11-1.15026429672983E-05*x12-3.54867343042830E-02*x13+
 5.74039763584371E-05*x14+ 3.58981136373263E-02*x15-1.60273110259053E-04*x16+
-3.64657746966166E-02*x17+ 3.42136847754702E-04*x18+ 4.84484213222674E-03*x21+
 1.04961436965878E-01*x22-1.93333579570378E-02*x23-1.73270432779827E-01*x24+
 4.33278434742597E-02*x25+ 2.39095330063948E-01*x26-7.65998823268819E-02*x27+
-3.01470084151145E-01*x28+ 3.51540715209933E-02;
 1.09666148073818E-01*x11*x22-6.87000856696835E-02*x11*x23+
-1.46268563004279E-01*x11*x24+ 1.67732815710094E-01*x11*x25+ 1.34582096848949E-01*x11*x26+
-2.76370603695369E-01*x11*x27-6.72919309015232E-02*x11*x28-1.09666148073818E-01*x12*x21+
 1.26087828516216E-01*x12*x23-1.17209231190407E-02*x12*x24-1.47857179787872E-01*x12*x25+
 3.85187670459295E-02*x12*x26+ 1.67571882683498E-01*x12*x27-8.18084558448740E-02*x12*x28+
 6.87000856696835E-02*x13*x21-1.26087828516216E-01*x13*x22+ 1.75513722783869E-01*x13*x24+
-1.00224780194645E-01*x13*x25-1.78864738930975E-01*x13*x26+ 2.12780032848162E-01*x13*x27+
 1.28617094854006E-01*x13*x28+ 1.46268563004279E-01*x14*x21+ 1.17209231190407E-02*x14*x22
-1.75513722783869E-01*x14*x23+ 2.15133302931685E-01*x14*x25-3.69909800370509E-02*x14*x26+
-2.53039042272502E-01*x14*x27+ 1.01920984056781E-01*x14*x28-1.67732815710094E-01*x15*x21+
 1.47857179787872E-01*x15*x22+ 1.00224780194645E-01*x15*x23-2.15133302931685E-01*x15*x24+
 2.40363968346693E-01*x15*x26-1.16317337176263E-01*x15*x27-2.15851091600313E-01*x15*x28+
-1.34582096848949E-01*x16*x21-3.85187670459295E-02*x16*x22+ 1.78864738930975E-01*x16*x23+
 3.69909800370509E-02*x16*x24-2.40363968346693E-01*x16*x25+ 3.02715385099606E-01*x16*x27+
-7.67597974850454E-02*x16*x28+ 2.76370603695369E-01*x17*x21-1.67571882683498E-01*x17*x22+
-2.12780032848162E-01*x17*x23+ 2.53039042272502E-01*x17*x24+ 1.16317337176263E-01*x17*x25+
-3.02715385099606E-01*x17*x26+ 3.08989496530880E-01*x17*x28+ 6.72919309015232E-02*x18*x21+
 8.18084558448740E-02*x18*x22-1.28617094854006E-01*x18*x23-1.01920984056781E-01*x18*x24+
 2.15851091600313E-01*x18*x25+ 7.67597974850454E-02*x18*x26-3.08989496530880E-01*x18*x27+
 4.03989735223750E-02*x11-1.24688987788719E-03*x12-4.56673047449227E-02*x13+
 5.98082512176471E-03*x14+ 5.25607404152185E-02*x15-1.57197798462696E-02*x16+
-5.85880734956333E-02*x17+ 3.09018172755273E-02*x18+ 2.45202776824656E-02*x21+
 1.03602528812701E-01*x22-9.30935512183135E-02*x23-1.35560453107760E-01*x24+
 1.91518014875911E-01*x25+ 1.18528414928598E-01*x26-2.98557057967123E-01*x27+
-4.52796805870746E-02*x28+ 3.84440411052168E-02;
 8.70610839178318E-02*x11*x22-9.28144441975650E-02*x11*x23+
-5.72885551402084E-02*x11*x24+ 1.81793241390571E-01*x11*x25-5.51837482137568E-02*x11*x26+
-2.01637769897399E-01*x11*x27+ 2.08096112924016E-01*x11*x28-8.70610839178318E-02*x12*x21+
 1.26375809016037E-01*x12*x23-5.09781480678748E-02*x12*x24-1.48589138848847E-01*x12*x25+
 1.43730448473742E-01*x12*x26+ 1.13359218630154E-01*x12*x27-2.42916913273311E-01*x12*x28+
 9.28144441975650E-02*x13*x21-1.26375809016037E-01*x13*x22+ 1.37505707988236E-01*x13*x24+
-1.05478236717080E-01*x13*x25-7.31253342804101E-02*x13*x26+ 1.71842145237536E-01*x13*x27+
-4.30975145652440E-02*x13*x28+ 5.72885551402084E-02*x14*x21+ 5.09781480678748E-02*x14*x22+
-1.37505707988236E-01*x14*x23+ 2.04223736386126E-01*x14*x25-1.26891080525973E-01*x14*x26+
-1.92661349735559E-01*x14*x27+ 2.81695475566136E-01*x14*x28-1.81793241390571E-01*x15*x21+
 1.48589138848847E-01*x15*x22+ 1.05478236717080E-01*x15*x23-2.04223736386126E-01*x15*x24+
 2.05941824776155E-01*x15*x26-1.07433107552462E-01*x15*x27-1.52075189505694E-01*x15*x28+
 5.51837482137568E-02*x16*x21-1.43730448473742E-01*x16*x22+ 7.31253342804101E-02*x16*x23+
 1.26891080525973E-01*x16*x24-2.05941824776155E-01*x16*x25+ 2.61033971727854E-01*x16*x27+
-1.89575882984077E-01*x16*x28+ 2.01637769897399E-01*x17*x21-1.13359218630154E-01*x17*x22+
-1.71842145237536E-01*x17*x23+ 1.92661349735559E-01*x17*x24+ 1.07433107552462E-01*x17*x25+
-2.61033971727854E-01*x17*x26+ 2.91652834528411E-01*x17*x28-2.08096112924016E-01*x18*x21+
 2.42916913273311E-01*x18*x22+ 4.30975145652440E-02*x18*x23-2.81695475566136E-01*x18*x24+
 1.52075189505694E-01*x18*x25+ 1.89575882984077E-01*x18*x26-2.91652834528411E-01*x18*x27+
 3.94056625817443E-02*x11-6.24517045064180E-03*x12-5.05424498095904E-02*x13+
 2.71832647496104E-02*x14+ 5.42139320707918E-02*x15-6.10968920005529E-02*x16+
-3.68446237196054E-02*x17+ 9.50218611145275E-02*x18+ 3.85912774161110E-02*x21+
 6.88184331926186E-02*x22-1.29384427861877E-01*x23-2.26874817253060E-02*x24+
 2.09565168440392E-01*x25-1.07331545616617E-01*x26-2.09635254114777E-01*x27+
 2.72168934359516E-01*x28+ 3.39169342981205E-02;
 7.47781044480180E-02*x11*x22-1.03965931164657E-01*x11*x23+
 1.55735916183540E-02*x11*x24+ 1.29673027259837E-01*x11*x25-1.58472645307897E-01*x11*x26+
-8.29778540500510E-03*x11*x27+ 2.08854787685352E-01*x11*x28-7.47781044480180E-02*x12*x21+
 1.36091194874051E-01*x12*x23-1.18926593447787E-01*x12*x24-9.38885560999291E-02*x12*x25+
 2.48486506175142E-01*x12*x26-9.75468195926671E-02*x12*x27-2.31148127842420E-01*x12*x28+
 1.03965931164657E-01*x13*x21-1.36091194874051E-01*x13*x22+ 1.37003811634985E-01*x13*x24+
-1.05460630759561E-01*x13*x25-5.70674982075000E-02*x13*x26+ 1.50723283839228E-01*x13*x27+
-5.87306577985588E-02*x13*x28-1.55735916183540E-02*x14*x21+ 1.18926593447787E-01*x14*x22+
-1.37003811634985E-01*x14*x23+ 1.86677496932823E-01*x14*x25-2.00283018696944E-01*x14*x26+
-3.35122386677441E-02*x14*x27+ 2.84021399513796E-01*x14*x28-1.29673027259837E-01*x15*x21+
 9.38885560999291E-02*x15*x22+ 1.05460630759561E-01*x15*x23-1.86677496932823E-01*x15*x24+
 2.31929249446437E-01*x15*x26-1.79574737629087E-01*x15*x27-1.38604784233871E-01*x15*x28+
 1.58472645307897E-01*x16*x21-2.48486506175142E-01*x16*x22+ 5.70674982075000E-02*x16*x23+
 2.00283018696944E-01*x16*x24-2.31929249446437E-01*x16*x25+ 2.34298400258090E-01*x16*x27+
-2.04163255078774E-01*x16*x28+ 8.29778540500510E-03*x17*x21+ 9.75468195926671E-02*x17*x22+
-1.50723283839228E-01*x17*x23+ 3.35122386677441E-02*x17*x24+ 1.79574737629087E-01*x17*x25+
-2.34298400258090E-01*x17*x26+ 2.98097123771152E-01*x17*x28-2.08854787685352E-01*x18*x21+
 2.31148127842420E-01*x18*x22+ 5.87306577985588E-02*x18*x23-2.84021399513796E-01*x18*x24+
 1.38604784233871E-01*x18*x25+ 2.04163255078774E-01*x18*x26-2.98097123771152E-01*x18*x27+
 4.49996439004562E-02*x11-1.85069481540755E-02*x12-5.61656816586600E-02*x13+
 6.77127982872445E-02*x14+ 2.44068716613928E-02*x15-1.10312494844420E-01*x16+
 6.07549343903777E-02*x17+ 8.74095266521031E-02*x18+ 5.51942508464642E-02*x21+
 3.33010653002596E-02*x22-1.46749210769555E-01*x23+ 9.47159796948122E-02*x24+
 1.27047329369764E-01*x25-2.53982560853732E-01*x26+ 6.83047340388718E-02*x27+
 2.63621748897238E-01*x28+ 3.36998809669997E-02;
 6.17544191308321E-02*x11*x22-9.31206181437554E-02*x11*x23+
 5.50538058077333E-02*x11*x24+ 4.29230152455850E-02*x11*x25-1.27451317851331E-01*x11*x26+
 1.13004878524310E-01*x11*x27+ 1.21249838675115E-02*x11*x28-6.17544191308321E-02*x12*x21+
 1.24530659024292E-01*x12*x23-1.73103189411473E-01*x12*x24+ 4.04716375488550E-02*x12*x25+
 1.91243383463927E-01*x12*x26-2.90495864846351E-01*x12*x27+ 1.06614414106906E-01*x12*x28+
 9.31206181437554E-02*x13*x21-1.24530659024292E-01*x13*x22+ 1.50006904963251E-01*x13*x24+
-1.47584179560868E-01*x13*x25-3.13678843682341E-02*x13*x26+ 2.10164433408536E-01*x13*x27+
-1.85216419121638E-01*x13*x28-5.50538058077333E-02*x14*x21+ 1.73103189411473E-01*x14*x22+
-1.50006904963251E-01*x14*x23+ 1.56397366349474E-01*x14*x25-1.86764828847131E-01*x14*x26+
 5.77869893465161E-02*x14*x27+ 1.29033723282119E-01*x14*x28-4.29230152455850E-02*x15*x21+
-4.04716375488550E-02*x15*x22+ 1.47584179560868E-01*x15*x23-1.56397366349474E-01*x15*x24+
 2.16452626279642E-01*x15*x26-2.75971358168906E-01*x15*x27+ 6.61571143772355E-02*x15*x28+
 1.27451317851331E-01*x16*x21-1.91243383463927E-01*x16*x22+ 3.13678843682341E-02*x16*x23+
 1.86764828847131E-01*x16*x24-2.16452626279642E-01*x16*x25+ 2.49579636648467E-01*x16*x27+
-2.57584327454215E-01*x16*x28-1.13004878524310E-01*x17*x21+ 2.90495864846351E-01*x17*x22+
-2.10164433408536E-01*x17*x23-5.77869893465161E-02*x17*x24+ 2.75971358168906E-01*x17*x25+
-2.49579636648467E-01*x17*x26+ 2.52131050847477E-01*x17*x28-1.21249838675115E-02*x18*x21+
-1.06614414106906E-01*x18*x22+ 1.85216419121638E-01*x18*x23-1.29033723282119E-01*x18*x24+
-6.61571143772355E-02*x18*x25+ 2.57584327454215E-01*x18*x26-2.52131050847477E-01*x18*x27+
 5.00462546184607E-02*x11-3.88990525396262E-02*x12-4.22640077968923E-02*x13+
 1.05605809263888E-01*x14-5.98356936199952E-02*x15-7.47036344688870E-02*x16+
 1.64238405178285E-01*x17-9.40386547332492E-02*x18+ 6.80907005283795E-02*x21+
-1.36224047417491E-02*x22-1.16766591955013E-01*x23+ 1.78720039161166E-01*x24+
-5.40925958928831E-02*x25-1.82751334619899E-01*x26+ 2.95374306837713E-01*x27+
-1.20228182617388E-01*x28+ 3.18505692187795E-02;
 5.06220113999399E-02*x11*x22-7.25233794493488E-02*x11*x23+
 5.17212948797265E-02*x11*x24+ 1.23968548196029E-03*x11*x25-5.60014680182141E-02*x11*x26+
 7.91778657330699E-02*x11*x27-5.44294920433267E-02*x11*x28-5.06220113999399E-02*x12*x21+
 8.73684801190553E-02*x12*x23-1.57938716118604E-01*x12*x24+ 1.40519311516436E-01*x12*x25+
-1.05113688851838E-02*x12*x26-1.68502712709901E-01*x12*x27+ 2.71061293822817E-01*x12*x28+
 7.25233794493488E-02*x13*x21-8.73684801190553E-02*x13*x22+ 1.37004404281524E-01*x13*x24+
-2.03453883017072E-01*x13*x25+ 1.11711940774171E-01*x13*x26+ 1.04751593966534E-01*x13*x27+
-2.94394842487108E-01*x13*x28-5.17212948797265E-02*x14*x21+ 1.57938716118604E-01*x14*x22+
-1.37004404281524E-01*x14*x23+ 1.47438532651002E-01*x14*x25-1.85462041306109E-01*x14*x26+
 7.48700389805950E-02*x14*x27+ 1.07129623379322E-01*x14*x28-1.23968548196029E-03*x15*x21+
-1.40519311516436E-01*x15*x22+ 2.03453883017072E-01*x15*x23-1.47438532651002E-01*x15*x24+
 1.55194484003375E-01*x15*x26-2.23912666314611E-01*x15*x27+ 1.57726358119397E-01*x15*x28+
 5.60014680182141E-02*x16*x21+ 1.05113688851838E-02*x16*x22-1.11711940774171E-01*x16*x23+
 1.85462041306109E-01*x16*x24-1.55194484003375E-01*x16*x25+ 2.02849842333208E-01*x16*x27+
-3.11168173892701E-01*x16*x28-7.91778657330699E-02*x17*x21+ 1.68502712709901E-01*x17*x22+
-1.04751593966534E-01*x17*x23-7.48700389805950E-02*x17*x24+ 2.23912666314611E-01*x17*x25+
-2.02849842333208E-01*x17*x26+ 2.42790385587628E-01*x17*x28+ 5.44294920433267E-02*x18*x21+
-2.71061293822817E-01*x18*x22+ 2.94394842487108E-01*x18*x23-1.07129623379322E-01*x18*x24+
-1.57726358119397E-01*x18*x25+ 3.11168173892701E-01*x18*x26-2.42790385587628E-01*x18*x27+
 4.97198948763377E-02*x11-6.10082421743452E-02*x12+ 1.59164457307720E-03*x13+
 9.27910793906784E-02*x14-1.39509202291483E-01*x15+ 7.78154634716169E-02*x16+
 7.00769222067556E-02*x17-2.00633896619934E-01*x18+ 7.01437608020686E-02*x21+
-5.89950508529677E-02*x22-3.65420745108934E-02*x23+ 1.58569659347986E-01*x24+
-1.96153570522525E-01*x25+ 7.98292340941004E-02*x26+ 1.41209556115732E-01*x27+
-3.12160411421906E-01*x28+ 2.65915909338954E-02;
 5.52052740993377E-02*x11*x22-8.20073139803128E-02*x11*x23+
 6.37172755995201E-02*x11*x24-7.28370499538614E-03*x11*x25-5.92993617762620E-02*x11*x26+
 1.00640771161039E-01*x11*x27-9.20921112605783E-02*x11*x28-5.52052740993377E-02*x12*x21+
 6.36358943374222E-02*x12*x23-1.17670501442470E-01*x12*x24+ 1.40523900246031E-01*x12*x25+
-1.15632212849186E-01*x12*x26+ 4.14816229696984E-02*x12*x27+ 6.24605738486459E-02*x12*x28+
 8.20073139803128E-02*x13*x21-6.36358943374222E-02*x13*x22+ 1.01351474687473E-01*x13*x24+
-2.00351917605564E-01*x13*x25+ 2.40126606051782E-01*x13*x26-1.77630889782092E-01*x13*x27+
 1.33708233956712E-02*x13*x28-6.37172755995201E-02*x14*x21+ 1.17670501442470E-01*x14*x22+
-1.01351474687473E-01*x14*x23+ 1.46665748757270E-01*x14*x25-2.59858418308211E-01*x14*x26+
 2.62394241265543E-01*x14*x27-1.24203845099501E-01*x14*x28+ 7.28370499538614E-03*x15*x21+
-1.40523900246031E-01*x15*x22+ 2.00351917605564E-01*x15*x23-1.46665748757270E-01*x15*x24+
 1.66201665963071E-01*x15*x26-2.61652058203112E-01*x15*x27+ 2.26177633680504E-01*x15*x28+
 5.92993617762620E-02*x16*x21+ 1.15632212849186E-01*x16*x22-2.40126606051782E-01*x16*x23+
 2.59858418308211E-01*x16*x24-1.66201665963071E-01*x16*x25+ 1.66242835569273E-01*x16*x27+
-2.59987601014195E-01*x16*x28-1.00640771161039E-01*x17*x21-4.14816229696984E-02*x17*x22+
 1.77630889782092E-01*x17*x23-2.62394241265543E-01*x17*x24+ 2.61652058203112E-01*x17*x25+
-1.66242835569273E-01*x17*x26+ 1.83066033489812E-01*x17*x28+ 9.20921112605783E-02*x18*x21+
-6.24605738486459E-02*x18*x22-1.33708233956712E-02*x18*x23+ 1.24203845099501E-01*x18*x24+
-2.26177633680504E-01*x18*x25+ 2.59987601014195E-01*x18*x26-1.83066033489812E-01*x18*x27+
 5.26902219386902E-02*x11-9.12167547115175E-02*x12+ 7.47655312889082E-02*x13+
 7.02834545654022E-03*x14-1.22086877923488E-01*x15+ 2.08345714854775E-01*x16+
-2.05882511101622E-01*x17+ 9.25506142044523E-02*x18+ 7.41155395919531E-02*x21+
-1.03560491210608E-01*x22+ 6.84048604782578E-02*x23+ 3.84495935108830E-02*x24+
-1.74995972445700E-01*x25+ 2.66482054904437E-01*x26-2.44484621931488E-01*x27+
 8.89010194192818E-02*x28+ 2.36200933801704E-02;
 6.21070804821762E-02*x11*x22-1.19460839651753E-01*x11*x23+
 1.51847459814976E-01*x11*x24-1.44178180839228E-01*x11*x25+ 9.08920979568109E-02*x11*x26+
 1.77826838177085E-03*x11*x27-1.16316815967881E-01*x11*x28-6.21070804821762E-02*x12*x21+
 6.48913422408734E-02*x12*x23-1.16602013146131E-01*x12*x24+ 1.41812838958244E-01*x12*x25+
-1.32978575516525E-01*x12*x26+ 9.02338991485880E-02*x12*x27-2.16343406290137E-02*x12*x28+
 1.19460839651753E-01*x13*x21-6.48913422408734E-02*x13*x22+ 6.56251892798338E-02*x13*x24+
-1.22130441166408E-01*x13*x25+ 1.60812937492014E-01*x13*x26-1.75419799077919E-01*x13*x27+
 1.63144213695607E-01*x13*x28-1.51847459814976E-01*x14*x21+ 1.16602013146131E-01*x14*x22+
-6.56251892798338E-02*x14*x23+ 7.60372761156268E-02*x14*x25-1.54479283628736E-01*x14*x26+
 2.23954144038926E-01*x14*x27-2.71271720439761E-01*x14*x28+ 1.44178180839228E-01*x15*x21+
-1.41812838958244E-01*x15*x22+ 1.22130441166408E-01*x15*x23-7.60372761156268E-02*x15*x24+
 1.01163065622469E-01*x15*x26-2.13533474991284E-01*x15*x27+ 3.15816129405170E-01*x15*x28+
-9.08920979568109E-02*x16*x21+ 1.32978575516525E-01*x16*x22-1.60812937492014E-01*x16*x23+
 1.54479283628736E-01*x16*x24-1.01163065622469E-01*x16*x25+ 1.35862448069089E-01*x16*x27+
-2.80709300265905E-01*x16*x28-1.77826838177085E-03*x17*x21-9.02338991485880E-02*x17*x22+
 1.75419799077919E-01*x17*x23-2.23954144038926E-01*x17*x24+ 2.13533474991284E-01*x17*x25+
-1.35862448069089E-01*x17*x26+ 1.68374492831468E-01*x17*x28+ 1.16316815967881E-01*x18*x21+
 2.16343406290137E-02*x18*x22-1.63144213695607E-01*x18*x23+ 2.71271720439761E-01*x18*x24+
-3.15816129405170E-01*x18*x25+ 2.80709300265905E-01*x18*x26-1.68374492831468E-01*x18*x27+
 4.69617571986567E-02*x11-1.06473377322043E-01*x12+ 1.55730836512497E-01*x13+
-1.72152295213730E-01*x14+ 1.39941495749199E-01*x15-5.52703659335981E-02*x16+
-7.12781967680265E-02*x17+ 2.15766556456721E-01*x18+ 6.43649193426506E-02*x21+
-1.22531028635551E-01*x22+ 1.68433831917568E-01*x23-1.78738819942904E-01*x24+
 1.37480763826661E-01*x25-4.15081010096656E-02*x26-9.70226043133603E-02*x27+
 2.51902223970491E-01*x28+ 1.76932793951543E-02;
 3.80763166686081E-02*x11*x22+ 8.41607820904365E-02*x11*x23+
 1.31920245452141E-01*x11*x24+ 1.74391082398050E-01*x11*x25+ 2.04712122107732E-01*x11*x26+
 2.16868567050131E-01*x11*x27+ 2.06376959977465E-01*x11*x28-3.80763166686081E-02*x12*x21+
 4.67337055845740E-02*x12*x23+ 9.63333232550634E-02*x12*x24+ 1.42390492988520E-01*x12*x25+
 1.78527635766573E-01*x12*x26+ 1.99058196602653E-01*x12*x27+ 1.99598756778582E-01*x12*x28+
-8.41607820904365E-02*x13*x21-4.67337055845740E-02*x13*x22+ 5.10124425063750E-02*x13*x24+
 1.00686045453334E-01*x13*x25+ 1.43345519846335E-01*x13*x26+ 1.73804147132172E-01*x13*x27+
 1.87876034596483E-01*x13*x28-1.31920245452141E-01*x14*x21-9.63333232550634E-02*x14*x22+
-5.10124425063750E-02*x14*x23+ 5.21194392019872E-02*x14*x25+ 1.00608746650142E-01*x14*x26+
 1.40983604803497E-01*x14*x27+ 1.69400276933406E-01*x14*x28-1.74391082398050E-01*x15*x21+
-1.42390492988520E-01*x15*x22-1.00686045453334E-01*x15*x23-5.21194392019872E-02*x15*x24+
 5.21207885973820E-02*x15*x26+ 1.00691257064772E-01*x15*x27+ 1.42401540948656E-01*x15*x28+
-2.04712122107732E-01*x16*x21-1.78527635766573E-01*x16*x22-1.43345519846335E-01*x16*x23+
-1.00608746650142E-01*x16*x24-5.21207885973820E-02*x16*x25+ 5.33820883887795E-02*x16*x27+
 1.05480116779528E-01*x16*x28-2.16868567050131E-01*x17*x21-1.99058196602653E-01*x17*x22+
-1.73804147132172E-01*x17*x23-1.40983604803497E-01*x17*x24-1.00691257064772E-01*x17*x25+
-5.33820883887795E-02*x17*x26+ 5.79276328985125E-02*x17*x28-2.06376959977465E-01*x18*x21+
-1.99598756778582E-01*x18*x22-1.87876034596483E-01*x18*x23-1.69400276933406E-01*x18*x24+
-1.42401540948656E-01*x18*x25-1.05480116779528E-01*x18*x26-5.79276328985125E-02*x18*x27+
 2.50194810815961E-02*x11+ 6.28969062735883E-02*x12+ 1.08314041933266E-01*x13+
 1.54614891096816E-01*x14+ 1.94507346481870E-01*x15+ 2.20848314317292E-01*x16+
 2.27438731187529E-01*x17+ 2.09752825072583E-01*x18-3.37328986508037E-02*x21+
-7.15094984012535E-02*x22-1.16656030491138E-01*x23-1.62409310885184E-01*x24+
-2.01368604685386E-01*x25-2.26298320814277E-01*x26-2.30940207627648E-01*x27+
-2.10757472319843E-01*x28+ 8.73415423203589E-03;
 5.02517061451570E-02*x11*x22+ 8.22742345300658E-02*x11*x23+
 7.81293687366289E-02*x11*x24+ 3.45795688546645E-02*x11*x25-3.44255250219385E-02*x11*x26+
-1.02217281170695E-01*x11*x27-1.39399785079317E-01*x11*x28-5.02517061451570E-02*x12*x21+
 5.10316485810222E-02*x12*x23+ 8.93034872554681E-02*x12*x24+ 1.04224133066684E-01*x12*x25+
 9.10135892286268E-02*x12*x26+ 5.22345394034698E-02*x12*x27-2.49862780181812E-03*x12*x28+
-8.22742345300658E-02*x13*x21-5.10316485810222E-02*x13*x22+ 6.68694821202251E-02*x13*x24+
 1.35523923149483E-01*x13*x25+ 1.83971160184117E-01*x13*x26+ 1.89324380143029E-01*x13*x27+
 1.37472509584302E-01*x13*x28-7.81293687366289E-02*x14*x21-8.93034872554681E-02*x14*x22+
-6.68694821202251E-02*x14*x23+ 1.00591403255221E-01*x14*x25+ 2.02682744313718E-01*x14*x26+
 2.62864930754803E-01*x14*x27+ 2.43845864297344E-01*x14*x28-3.45795688546645E-02*x15*x21+
-1.04224133066684E-01*x15*x22-1.35523923149483E-01*x15*x23-1.00591403255221E-01*x15*x24+
 1.34028905541301E-01*x15*x26+ 2.47946912097908E-01*x15*x27+ 2.87401590619440E-01*x15*x28+
 3.44255250219385E-02*x16*x21-9.10135892286268E-02*x16*x22-1.83971160184117E-01*x16*x23+
-2.02682744313718E-01*x16*x24-1.34028905541301E-01*x16*x25+ 1.49347370925362E-01*x16*x27+
 2.54186222350003E-01*x16*x28+ 1.02217281170695E-01*x17*x21-5.22345394034698E-02*x17*x22+
-1.89324380143029E-01*x17*x23-2.62864930754803E-01*x17*x24-2.47946912097908E-01*x17*x25+
-1.49347370925362E-01*x17*x26+ 1.49982698803284E-01*x17*x28+ 1.39399785079317E-01*x18*x21+
 2.49862780181812E-03*x18*x22-1.37472509584302E-01*x18*x23-2.43845864297344E-01*x18*x24+
-2.87401590619440E-01*x18*x25-2.54186222350003E-01*x18*x26-1.49982698803284E-01*x18*x27+
 4.40265179181883E-02*x11+ 8.66811783551952E-02*x12+ 9.72082777165868E-02*x13+
 5.65279148179278E-02*x14-3.16651514683834E-02*x15-1.39120778728911E-01*x16+
-2.22082395230492E-01*x17-2.38267168817362E-01*x18-6.13440148846316E-02*x21+
-9.99514493232346E-02*x22-1.01348653882254E-01*x23-4.63846776018603E-02*x24+
 5.84507267812891E-02*x25+ 1.79576392211120E-01*x26+ 2.67076340110360E-01*x27+
 2.74218245509010E-01*x28+ 1.82452953539884E-02;
 5.56228200642702E-02*x11*x22+ 7.87028715085628E-02*x11*x23+
 5.56907228560392E-02*x11*x24+ 8.70602425792709E-06*x11*x25-5.57732009702447E-02*x11*x26+
-7.89741811139716E-02*x11*x27-5.59331698756202E-02*x11*x28-5.56228200642702E-02*x12*x21+
 8.11117125987734E-02*x12*x23+ 1.53453662341955E-01*x12*x24+ 1.66493200852263E-01*x12*x25+
 8.75066268871735E-02*x12*x26-6.54075567659762E-02*x12*x27-2.18356868962363E-01*x12*x28+
-7.87028715085628E-02*x13*x21-8.11117125987734E-02*x13*x22+ 1.35916768591827E-01*x13*x24+
 2.35564950105869E-01*x13*x25+ 2.05147503249138E-01*x13*x26+ 2.26160511808778E-02*x13*x27+
-2.27397269441793E-01*x13*x28-5.56907228560392E-02*x14*x21-1.53453662341955E-01*x14*x22+
-1.35916768591827E-01*x14*x23+ 1.66672432718431E-01*x14*x25+ 2.41481989585178E-01*x14*x26+
 1.52388591517962E-01*x14*x27-6.43135695952369E-02*x14*x28-8.70602425792709E-06*x15*x21+
-1.66493200852263E-01*x15*x22-2.35564950105869E-01*x15*x23-1.66672432718431E-01*x15*x24+
 1.66957025469019E-01*x15*x26+ 2.36379506529577E-01*x15*x27+ 1.67387979527994E-01*x15*x28+
 5.57732009702447E-02*x16*x21-8.75066268871735E-02*x16*x22-2.05147503249138E-01*x16*x23+
-2.41481989585178E-01*x16*x24-1.66957025469019E-01*x16*x25+ 1.89827718133341E-01*x16*x27+
 3.06942088572043E-01*x16*x28+ 7.89741811139716E-02*x17*x21+ 6.54075567659762E-02*x17*x22+
-2.26160511808778E-02*x17*x23-1.52388591517962E-01*x17*x24-2.36379506529577E-01*x17*x25+
-1.89827718133341E-01*x17*x26+ 2.44254119396940E-01*x17*x28+ 5.59331698756202E-02*x18*x21+
 2.18356868962363E-01*x18*x22+ 2.27397269441793E-01*x18*x23+ 6.43135695952369E-02*x18*x24+
-1.67387979527994E-01*x18*x25-3.06942088572043E-01*x18*x26-2.44254119396940E-01*x18*x27+
 5.55911499757188E-02*x11+ 8.04997971983394E-02*x12+ 3.28367352354669E-02*x13+
-7.27682209776133E-02*x14-1.66385804506762E-01*x15-1.68174238119399E-01*x16+
-4.89246007791860E-02*x17+ 1.37283593500727E-01*x18-7.86122259316262E-02*x21+
-8.63544149527502E-02*x22-7.55028506685289E-03*x23+ 1.30417590750527E-01*x24+
 2.35292804333259E-01*x25+ 2.10261774747902E-01*x26+ 3.01663161703806E-02*x27+
-2.21769470538829E-01*x28+ 2.74658316557292E-02;
 5.30120656385646E-02*x11*x22+ 7.87781756792552E-02*x11*x23+
 5.39662221740742E-02*x11*x24-1.47478173763072E-02*x11*x25-8.49747320304881E-02*x11*x26+
-1.01957442584773E-01*x11*x27-4.11581250206766E-02*x11*x28-5.30120656385646E-02*x12*x21+
 1.03664180443218E-01*x12*x23+ 1.66095549388702E-01*x12*x24+ 9.34226101199596E-02*x12*x25+
-9.83355259538747E-02*x12*x26-2.60770213615841E-01*x12*x27-2.27554332040879E-01*x12*x28+
-7.87781756792552E-02*x13*x21-1.03664180443218E-01*x13*x22+ 1.41295006820949E-01*x13*x24+
 1.67669059619133E-01*x13*x25+ 2.00358654661370E-02*x13*x26-1.88139565104084E-01*x13*x27+
-2.57671375048381E-01*x13*x28-5.39662221740742E-02*x14*x21-1.66095549388702E-01*x14*x22+
-1.41295006820949E-01*x14*x23+ 1.41311455665444E-01*x14*x25+ 1.66134404563242E-01*x14*x26+
 5.39857129040022E-02*x14*x27-1.02694852330232E-01*x14*x28+ 1.47478173763072E-02*x15*x21+
-9.34226101199596E-02*x15*x22-1.67669059619133E-01*x15*x23-1.41311455665444E-01*x15*x24+
 1.77106768540502E-01*x15*x26+ 2.52224125469382E-01*x15*x27+ 1.35837551554738E-01*x15*x28+
 8.49747320304881E-02*x16*x21+ 9.83355259538747E-02*x16*x22-2.00358654661370E-02*x16*x23+
-1.66134404563242E-01*x16*x24-1.77106768540502E-01*x16*x25+ 2.28869411819851E-01*x16*x27+
 2.88407220737396E-01*x16*x28+ 1.01957442584773E-01*x17*x21+ 2.60770213615841E-01*x17*x22+
 1.88139565104084E-01*x17*x23-5.39857129040022E-02*x17*x24-2.52224125469382E-01*x17*x25+
-2.28869411819851E-01*x17*x26+ 2.35192583803692E-01*x17*x28+ 4.11581250206766E-02*x18*x21+
 2.27554332040879E-01*x18*x22+ 2.57671375048381E-01*x18*x23+ 1.02694852330232E-01*x18*x24+
-1.35837551554738E-01*x18*x25-2.88407220737396E-01*x18*x26-2.35192583803692E-01*x18*x27+
 4.78525115318335E-02*x11+ 4.64460042733650E-02*x12-2.45540310525770E-02*x13+
-1.02647835756774E-01*x14-9.72511381113145E-02*x15+ 1.43149510123702E-02*x16+
 1.46060972100363E-01*x17+ 1.69346652302113E-01*x18-6.65903666558965E-02*x21+
-3.34033978145379E-02*x22+ 8.05774495301260E-02*x23+ 1.74633986326751E-01*x24+
 1.26644245827537E-01*x25-6.99794264887512E-02*x26-2.63318528615275E-01*x27+
-2.59905080388915E-01*x28+ 2.81901857084959E-02;
 7.19386105014067E-02*x11*x22+ 1.05983357143432E-01*x11*x23+
 4.08051528126271E-02*x11*x24-9.38974984819799E-02*x11*x25-1.66418093164157E-01*x11*x26+
-7.06999395977262E-02*x11*x27+ 1.29793553565229E-01*x11*x28-7.19386105014067E-02*x12*x21+
 1.40302546263150E-01*x12*x23+ 1.55285605654952E-01*x12*x24-4.11289567082487E-02*x12*x25+
-2.60513786025372E-01*x12*x26-2.18342752576973E-01*x12*x27+ 1.08405743242592E-01*x12*x28+
-1.05983357143432E-01*x13*x21-1.40302546263150E-01*x13*x22+ 1.49191413182556E-01*x13*x24+
 1.22536050607242E-01*x13*x25-5.92344414247729E-02*x13*x26-1.83786096057204E-01*x13*x27+
-9.34291252558175E-02*x13*x28-4.08051528126271E-02*x14*x21-1.55285605654952E-01*x14*x22+
-1.49191413182556E-01*x14*x23+ 1.79356488370079E-01*x14*x25+ 2.11458484315311E-01*x14*x26+
 2.87630458990859E-02*x14*x27-2.18680310156908E-01*x14*x28+ 9.38974984819799E-02*x15*x21+
 4.11289567082487E-02*x15*x22-1.22536050607242E-01*x15*x23-1.79356488370079E-01*x15*x24+
 2.44889220904396E-01*x15*x26+ 2.44569965989445E-01*x15*x27-6.72900773621488E-02*x15*x28+
 1.66418093164157E-01*x16*x21+ 2.60513786025372E-01*x16*x22+ 5.92344414247729E-02*x16*x23+
-2.11458484315311E-01*x16*x24-2.44889220904396E-01*x16*x25+ 2.49071749910349E-01*x16*x27+
 2.19247117126079E-01*x16*x28+ 7.06999395977262E-02*x17*x21+ 2.18342752576973E-01*x17*x22+
 1.83786096057204E-01*x17*x23-2.87630458990859E-02*x17*x24-2.44569965989445E-01*x17*x25+
-2.49071749910349E-01*x17*x26+ 2.87400633801186E-01*x17*x28-1.29793553565229E-01*x18*x21+
-1.08405743242592E-01*x18*x22+ 9.34291252558175E-02*x18*x23+ 2.18680310156908E-01*x18*x24+
 6.72900773621488E-02*x18*x25-2.19247117126079E-01*x18*x26-2.87400633801186E-01*x18*x27+
 4.96878217664269E-02*x11+ 2.79794728369102E-02*x12-5.56859441935113E-02*x13+
-9.13848180570500E-02*x14-8.11239241075434E-03*x15+ 1.15210343829295E-01*x16+
 1.23311093620809E-01*x17-2.43942721450851E-02*x18-6.45483251424406E-02*x21+
 1.38376001193583E-02*x22+ 1.46275409222122E-01*x23+ 1.47182035869042E-01*x24+
-5.49652165786240E-02*x25-2.65762091553211E-01*x26-2.09511642960567E-01*x27+
 1.22235339175268E-01*x28+ 3.46627256353244E-02]

L=1e4 # Squared radius of a ball containing at least one real root
k=1 # relaxed order



sol=ASC_PolySys(x,h,k,L,method="LMBM",EigAlg="Arpack",tol=1e-3);
end