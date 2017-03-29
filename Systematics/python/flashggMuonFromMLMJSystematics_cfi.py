import FWCore.ParameterSet.Config as cms

binID = cms.PSet(
		variables = cms.vstring("abs(eta)","pt"),
		bins = cms.VPSet(
			# HighPtID SFs for Moriond2017
			# uncertainties are: stat (+) syst. with syst = 1% 
			# Preliminary numbers merged for 2016BCDEF and GH taken from : https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
			cms.PSet(lowBounds = cms.vdouble(0.0,0.00), upBounds = cms.vdouble(0.9,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			cms.PSet(lowBounds = cms.vdouble(0.0,20.00), upBounds = cms.vdouble(0.9,25.00), values = cms.vdouble(0.986792679377), uncertainties = cms.vdouble(0.00747498847049,0.00747498847049)),
			cms.PSet(lowBounds = cms.vdouble(0.0,25.00), upBounds = cms.vdouble(0.9,30.00), values = cms.vdouble(0.983975321203), uncertainties = cms.vdouble(0.0284405704478,0.0284405704478)),
			cms.PSet(lowBounds = cms.vdouble(0.0,30.00), upBounds = cms.vdouble(0.9,40.00), values = cms.vdouble(0.985440559068), uncertainties = cms.vdouble(0.0071116498314,0.0071116498314)),
			cms.PSet(lowBounds = cms.vdouble(0.0,40.00), upBounds = cms.vdouble(0.9,50.00), values = cms.vdouble(0.986990914465), uncertainties = cms.vdouble(0.0199690369541,0.0199690369541)),
			cms.PSet(lowBounds = cms.vdouble(0.0,50.00), upBounds = cms.vdouble(0.9,55.00), values = cms.vdouble(0.983094644476), uncertainties = cms.vdouble(0.0071406592439,0.0071406592439)),
			cms.PSet(lowBounds = cms.vdouble(0.0,55.00), upBounds = cms.vdouble(0.9,60.00), values = cms.vdouble(0.982697185499), uncertainties = cms.vdouble(0.00720095336476,0.00720095336476)),
			cms.PSet(lowBounds = cms.vdouble(0.0,60.00), upBounds = cms.vdouble(0.9,120.00), values = cms.vdouble(0.993846213993), uncertainties = cms.vdouble(0.00722924074947,0.00722924074947)),
			cms.PSet(lowBounds = cms.vdouble(0.0,120.00), upBounds = cms.vdouble(0.9,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),

			cms.PSet(lowBounds = cms.vdouble(0.9,0.00), upBounds = cms.vdouble(1.2,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			cms.PSet(lowBounds = cms.vdouble(0.9,20.00), upBounds = cms.vdouble(1.2,25.00), values = cms.vdouble(0.97358262229), uncertainties = cms.vdouble(0.014149852681,0.014149852681)),
			cms.PSet(lowBounds = cms.vdouble(0.9,25.00), upBounds = cms.vdouble(1.2,30.00), values = cms.vdouble(0.972545424629), uncertainties = cms.vdouble(0.0184129524702,0.0184129524702)),
			cms.PSet(lowBounds = cms.vdouble(0.9,30.00), upBounds = cms.vdouble(1.2,40.00), values = cms.vdouble(0.975390818992), uncertainties = cms.vdouble(0.00726695178253,0.00726695178253)),
			cms.PSet(lowBounds = cms.vdouble(0.9,40.00), upBounds = cms.vdouble(1.2,50.00), values = cms.vdouble(0.97708593518), uncertainties = cms.vdouble(0.00711334159536,0.00711334159536)),
			cms.PSet(lowBounds = cms.vdouble(0.9,50.00), upBounds = cms.vdouble(1.2,55.00), values = cms.vdouble(0.976316712692), uncertainties = cms.vdouble(0.00719504272429,0.00719504272429)),
			cms.PSet(lowBounds = cms.vdouble(0.9,55.00), upBounds = cms.vdouble(1.2,60.00), values = cms.vdouble(0.976151167477), uncertainties = cms.vdouble(0.00795625739171,0.00795625739171)),
			cms.PSet(lowBounds = cms.vdouble(0.9,60.00), upBounds = cms.vdouble(1.2,120.00), values = cms.vdouble(0.977022227485), uncertainties = cms.vdouble(0.00738405770892,0.00738405770892)),
			cms.PSet(lowBounds = cms.vdouble(0.9,120.00), upBounds = cms.vdouble(1.2,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),

			cms.PSet(lowBounds = cms.vdouble(1.2,0.00), upBounds = cms.vdouble(2.1,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			cms.PSet(lowBounds = cms.vdouble(1.2,20.00), upBounds = cms.vdouble(2.1,25.00), values = cms.vdouble(0.986872821242), uncertainties = cms.vdouble(0.00721311004636,0.00721311004636)),
			cms.PSet(lowBounds = cms.vdouble(1.2,25.00), upBounds = cms.vdouble(2.1,30.00), values = cms.vdouble(0.986750531472), uncertainties = cms.vdouble(0.00719484280123,0.00719484280123)),
			cms.PSet(lowBounds = cms.vdouble(1.2,30.00), upBounds = cms.vdouble(2.1,40.00), values = cms.vdouble(0.988776470703), uncertainties = cms.vdouble(0.00711242990378,0.00711242990378)),
			cms.PSet(lowBounds = cms.vdouble(1.2,40.00), upBounds = cms.vdouble(2.1,50.00), values = cms.vdouble(0.990218614475), uncertainties = cms.vdouble(0.00710934816151,0.00710934816151)),
			cms.PSet(lowBounds = cms.vdouble(1.2,50.00), upBounds = cms.vdouble(2.1,55.00), values = cms.vdouble(0.986546873029), uncertainties = cms.vdouble(0.0111281952296,0.0111281952296)),
			cms.PSet(lowBounds = cms.vdouble(1.2,55.00), upBounds = cms.vdouble(2.1,60.00), values = cms.vdouble(0.986223325783), uncertainties = cms.vdouble(0.00721970513803,0.00721970513803)),
			cms.PSet(lowBounds = cms.vdouble(1.2,60.00), upBounds = cms.vdouble(2.1,120.00), values = cms.vdouble(0.990005820233), uncertainties = cms.vdouble(0.00786295546941,0.00786295546941)),
			cms.PSet(lowBounds = cms.vdouble(1.2,120.00), upBounds = cms.vdouble(2.1,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),

			cms.PSet(lowBounds = cms.vdouble(2.1,0.00), upBounds = cms.vdouble(2.4,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			cms.PSet(lowBounds = cms.vdouble(2.1,20.00), upBounds = cms.vdouble(2.4,25.00), values = cms.vdouble(0.977910767205), uncertainties = cms.vdouble(0.00804485443992,0.00804485443992)),
			cms.PSet(lowBounds = cms.vdouble(2.1,25.00), upBounds = cms.vdouble(2.4,30.00), values = cms.vdouble(0.974876109599), uncertainties = cms.vdouble(0.0205533733304,0.0205533733304)),
			cms.PSet(lowBounds = cms.vdouble(2.1,30.00), upBounds = cms.vdouble(2.4,40.00), values = cms.vdouble(0.97495211615), uncertainties = cms.vdouble(0.00715515860995,0.00715515860995)),
			cms.PSet(lowBounds = cms.vdouble(2.1,40.00), upBounds = cms.vdouble(2.4,50.00), values = cms.vdouble(0.975918132749), uncertainties = cms.vdouble(0.00712816359643,0.00712816359643)),
			cms.PSet(lowBounds = cms.vdouble(2.1,50.00), upBounds = cms.vdouble(2.4,55.00), values = cms.vdouble(0.970796501484), uncertainties = cms.vdouble(0.00741039722618,0.00741039722618)),
			cms.PSet(lowBounds = cms.vdouble(2.1,55.00), upBounds = cms.vdouble(2.4,60.00), values = cms.vdouble(0.966171651071), uncertainties = cms.vdouble(0.00798936866738,0.00798936866738)),
			cms.PSet(lowBounds = cms.vdouble(2.1,60.00), upBounds = cms.vdouble(2.4,120.00), values = cms.vdouble(0.967327162846), uncertainties = cms.vdouble(0.00812426767543,0.00812426767543)),
			cms.PSet(lowBounds = cms.vdouble(2.1,120.00), upBounds = cms.vdouble(2.4,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000))
		)
)	


binIso = cms.PSet(
		variables = cms.vstring("abs(eta)","pt"),
		bins = cms.VPSet(
			# tkLooseISO_highptID SFs for Moriond2017
			# uncertainties are: stat (+) syst. with syst = 1% 
			# Preliminary numbers merged for 2016BCDEF and GH taken from : https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
			cms.PSet(lowBounds = cms.vdouble(0.0,0.00), upBounds = cms.vdouble(0.9,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			cms.PSet(lowBounds = cms.vdouble(0.0,20.00), upBounds = cms.vdouble(0.9,25.00), values = cms.vdouble(0.990820979203), uncertainties = cms.vdouble(0.0071850680993,0.0071850680993)),
			cms.PSet(lowBounds = cms.vdouble(0.0,25.00), upBounds = cms.vdouble(0.9,30.00), values = cms.vdouble(0.994687947027), uncertainties = cms.vdouble(0.00712489310061,0.00712489310061)),
			cms.PSet(lowBounds = cms.vdouble(0.0,30.00), upBounds = cms.vdouble(0.9,40.00), values = cms.vdouble(0.996669339513), uncertainties = cms.vdouble(0.00710763561077,0.00710763561077)),
			cms.PSet(lowBounds = cms.vdouble(0.0,40.00), upBounds = cms.vdouble(0.9,50.00), values = cms.vdouble(0.997877878982), uncertainties = cms.vdouble(0.00710631633081,0.00710631633081)),
			cms.PSet(lowBounds = cms.vdouble(0.0,50.00), upBounds = cms.vdouble(0.9,55.00), values = cms.vdouble(0.997917093292), uncertainties = cms.vdouble(0.00710844642708,0.00710844642708)),
			cms.PSet(lowBounds = cms.vdouble(0.0,55.00), upBounds = cms.vdouble(0.9,60.00), values = cms.vdouble(0.997972855217), uncertainties = cms.vdouble(0.00711114419668,0.00711114419668)),
			cms.PSet(lowBounds = cms.vdouble(0.0,60.00), upBounds = cms.vdouble(0.9,120.00), values = cms.vdouble(0.998730942045), uncertainties = cms.vdouble(0.00710912356167,0.00710912356167)),
			cms.PSet(lowBounds = cms.vdouble(0.0,120.00), upBounds = cms.vdouble(0.9,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),

			cms.PSet(lowBounds = cms.vdouble(0.9,0.00), upBounds = cms.vdouble(1.2,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			cms.PSet(lowBounds = cms.vdouble(0.9,20.00), upBounds = cms.vdouble(1.2,25.00), values = cms.vdouble(0.994909170637), uncertainties = cms.vdouble(0.00729354765757,0.00729354765757)),
			cms.PSet(lowBounds = cms.vdouble(0.9,25.00), upBounds = cms.vdouble(1.2,30.00), values = cms.vdouble(0.998491509847), uncertainties = cms.vdouble(0.00717071447953,0.00717071447953)),
			cms.PSet(lowBounds = cms.vdouble(0.9,30.00), upBounds = cms.vdouble(1.2,40.00), values = cms.vdouble(1.00010005997), uncertainties = cms.vdouble(0.00711238191705,0.00711238191705)),
			cms.PSet(lowBounds = cms.vdouble(0.9,40.00), upBounds = cms.vdouble(1.2,50.00), values = cms.vdouble(0.999207915723), uncertainties = cms.vdouble(0.00830962058259,0.00830962058259)),
			cms.PSet(lowBounds = cms.vdouble(0.9,50.00), upBounds = cms.vdouble(1.2,55.00), values = cms.vdouble(0.999246060343), uncertainties = cms.vdouble(0.00836390014695,0.00836390014695)),
			cms.PSet(lowBounds = cms.vdouble(0.9,55.00), upBounds = cms.vdouble(1.2,60.00), values = cms.vdouble(0.998718769954), uncertainties = cms.vdouble(0.00711103691333,0.00711103691333)),
			cms.PSet(lowBounds = cms.vdouble(0.9,60.00), upBounds = cms.vdouble(1.2,120.00), values = cms.vdouble(0.999313937016), uncertainties = cms.vdouble(0.0071165794739,0.0071165794739)),
			cms.PSet(lowBounds = cms.vdouble(0.9,120.00), upBounds = cms.vdouble(1.2,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),

			cms.PSet(lowBounds = cms.vdouble(1.2,0.00), upBounds = cms.vdouble(2.1,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			cms.PSet(lowBounds = cms.vdouble(1.2,20.00), upBounds = cms.vdouble(2.1,25.00), values = cms.vdouble(0.997260352619), uncertainties = cms.vdouble(0.0071524646053,0.0071524646053)),
			cms.PSet(lowBounds = cms.vdouble(1.2,25.00), upBounds = cms.vdouble(2.1,30.00), values = cms.vdouble(1.00083998444), uncertainties = cms.vdouble(0.00712293237365,0.00712293237365)),
			cms.PSet(lowBounds = cms.vdouble(1.2,30.00), upBounds = cms.vdouble(2.1,40.00), values = cms.vdouble(1.00079772922), uncertainties = cms.vdouble(0.00710818204618,0.00710818204618)),
			cms.PSet(lowBounds = cms.vdouble(1.2,40.00), upBounds = cms.vdouble(2.1,50.00), values = cms.vdouble(1.00007372819), uncertainties = cms.vdouble(0.00710637897594,0.00710637897594)),
			cms.PSet(lowBounds = cms.vdouble(1.2,50.00), upBounds = cms.vdouble(2.1,55.00), values = cms.vdouble(0.999467448457), uncertainties = cms.vdouble(0.00710877337951,0.00710877337951)),
			cms.PSet(lowBounds = cms.vdouble(1.2,55.00), upBounds = cms.vdouble(2.1,60.00), values = cms.vdouble(1.00008129962), uncertainties = cms.vdouble(0.00711263871536,0.00711263871536)),
			cms.PSet(lowBounds = cms.vdouble(1.2,60.00), upBounds = cms.vdouble(2.1,120.00), values = cms.vdouble(0.999228284006), uncertainties = cms.vdouble(0.00710950699049,0.00710950699049)),
			cms.PSet(lowBounds = cms.vdouble(1.2,120.00), upBounds = cms.vdouble(2.1,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),

			cms.PSet(lowBounds = cms.vdouble(2.1,0.00), upBounds = cms.vdouble(2.4,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			cms.PSet(lowBounds = cms.vdouble(2.1,20.00), upBounds = cms.vdouble(2.4,25.00), values = cms.vdouble(1.00031542892), uncertainties = cms.vdouble(0.00721554294369,0.00721554294369)),
			cms.PSet(lowBounds = cms.vdouble(2.1,25.00), upBounds = cms.vdouble(2.4,30.00), values = cms.vdouble(1.00174457993), uncertainties = cms.vdouble(0.00714628187831,0.00714628187831)),
			cms.PSet(lowBounds = cms.vdouble(2.1,30.00), upBounds = cms.vdouble(2.4,40.00), values = cms.vdouble(1.0004314503), uncertainties = cms.vdouble(0.00711181549784,0.00711181549784)),
			cms.PSet(lowBounds = cms.vdouble(2.1,40.00), upBounds = cms.vdouble(2.4,50.00), values = cms.vdouble(1.00001432953), uncertainties = cms.vdouble(0.00710662671288,0.00710662671288)),
			cms.PSet(lowBounds = cms.vdouble(2.1,50.00), upBounds = cms.vdouble(2.4,55.00), values = cms.vdouble(0.999956981123), uncertainties = cms.vdouble(0.00713230338238,0.00713230338238)),
			cms.PSet(lowBounds = cms.vdouble(2.1,55.00), upBounds = cms.vdouble(2.4,60.00), values = cms.vdouble(1.00055772784), uncertainties = cms.vdouble(0.0071287284802,0.0071287284802)),
			cms.PSet(lowBounds = cms.vdouble(2.1,60.00), upBounds = cms.vdouble(2.4,120.00), values = cms.vdouble(1.00045987414), uncertainties = cms.vdouble(0.00712572283075,0.00712572283075)),
			cms.PSet(lowBounds = cms.vdouble(2.1,120.00), upBounds = cms.vdouble(2.4,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000))
		)
)


binTrigger = cms.PSet(
		variables = cms.vstring("abs(eta)","pt"),
		bins = cms.VPSet(
			# Mu50_OR_TkMu50 trigger SFs for Moriond2017
			# uncertainties are: stat (+) syst. with syst = 1% 
			# Preliminary numbers merged for 2016BCDEF and GH taken from : https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
			cms.PSet(lowBounds = cms.vdouble(0.0,0.00), upBounds = cms.vdouble(0.9,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			cms.PSet(lowBounds = cms.vdouble(0.0,52.00), upBounds = cms.vdouble(0.9,55.00), values = cms.vdouble(0.976918626141), uncertainties = cms.vdouble(0.00713301391063,0.00713301391063)),
			cms.PSet(lowBounds = cms.vdouble(0.0,55.00), upBounds = cms.vdouble(0.9,60.00), values = cms.vdouble(0.981657641158), uncertainties = cms.vdouble(0.00713569925415,0.00713569925415)),
			cms.PSet(lowBounds = cms.vdouble(0.0,60.00), upBounds = cms.vdouble(0.9,80.00), values = cms.vdouble(0.980945665477), uncertainties = cms.vdouble(0.00713214355006,0.00713214355006)),
			cms.PSet(lowBounds = cms.vdouble(0.0,80.00), upBounds = cms.vdouble(0.9,120.00), values = cms.vdouble(0.978871910597), uncertainties = cms.vdouble(0.00726935141737,0.00726935141737)),
			cms.PSet(lowBounds = cms.vdouble(0.0,120.00), upBounds = cms.vdouble(0.9,200.00), values = cms.vdouble(0.971546057023), uncertainties = cms.vdouble(0.00755507754235,0.00755507754235)),
			cms.PSet(lowBounds = cms.vdouble(0.0,200.00), upBounds = cms.vdouble(0.9,300.00), values = cms.vdouble(0.977208905196), uncertainties = cms.vdouble(0.00985642270769,0.00985642270769)),
			cms.PSet(lowBounds = cms.vdouble(0.0,300.00), upBounds = cms.vdouble(0.9,400.00), values = cms.vdouble(1.00520971939), uncertainties = cms.vdouble(0.0184154737777,0.0184154737777)),
			cms.PSet(lowBounds = cms.vdouble(0.0,400.00), upBounds = cms.vdouble(0.9,800.00), values = cms.vdouble(0.964794536058), uncertainties = cms.vdouble(0.0354231365919,0.0354231365919)),
			cms.PSet(lowBounds = cms.vdouble(0.0,800.00), upBounds = cms.vdouble(0.9,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),

			cms.PSet(lowBounds = cms.vdouble(0.9,0.00), upBounds = cms.vdouble(1.2,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),		
			cms.PSet(lowBounds = cms.vdouble(0.9,52.00), upBounds = cms.vdouble(1.2,55.00), values = cms.vdouble(0.952914384419), uncertainties = cms.vdouble(0.007149529456,0.007149529456)),
			cms.PSet(lowBounds = cms.vdouble(0.9,55.00), upBounds = cms.vdouble(1.2,60.00), values = cms.vdouble(0.958066480363), uncertainties = cms.vdouble(0.00717378231989,0.00717378231989)),
			cms.PSet(lowBounds = cms.vdouble(0.9,60.00), upBounds = cms.vdouble(1.2,80.00), values = cms.vdouble(0.95763756959), uncertainties = cms.vdouble(0.00716157339882,0.00716157339882)),
			cms.PSet(lowBounds = cms.vdouble(0.9,80.00), upBounds = cms.vdouble(1.2,120.00), values = cms.vdouble(0.953281720938), uncertainties = cms.vdouble(0.00747026796001,0.00747026796001)),
			cms.PSet(lowBounds = cms.vdouble(0.9,120.00), upBounds = cms.vdouble(1.2,200.00), values = cms.vdouble(0.940346424224), uncertainties = cms.vdouble(0.0103033975804,0.0103033975804)),
			cms.PSet(lowBounds = cms.vdouble(0.9,200.00), upBounds = cms.vdouble(1.2,300.00), values = cms.vdouble(0.943458228056), uncertainties = cms.vdouble(0.0199257906567,0.0199257906567)),
			cms.PSet(lowBounds = cms.vdouble(0.9,300.00), upBounds = cms.vdouble(1.2,400.00), values = cms.vdouble(0.942463210436), uncertainties = cms.vdouble(0.0395713515752,0.0395713515752)),
			cms.PSet(lowBounds = cms.vdouble(0.9,400.00), upBounds = cms.vdouble(1.2,800.00), values = cms.vdouble(0.975455933657), uncertainties = cms.vdouble(0.0583310146172,0.0583310146172)),
			cms.PSet(lowBounds = cms.vdouble(0.9,800.00), upBounds = cms.vdouble(1.2,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			
			cms.PSet(lowBounds = cms.vdouble(1.2,0.00), upBounds = cms.vdouble(2.1,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			cms.PSet(lowBounds = cms.vdouble(1.2,52.00), upBounds = cms.vdouble(2.1,55.00), values = cms.vdouble(0.986467439092), uncertainties = cms.vdouble(0.00718225989039,0.00718225989039)),
			cms.PSet(lowBounds = cms.vdouble(1.2,55.00), upBounds = cms.vdouble(2.1,60.00), values = cms.vdouble(0.992607911292), uncertainties = cms.vdouble(0.00719339972107,0.00719339972107)),
			cms.PSet(lowBounds = cms.vdouble(1.2,60.00), upBounds = cms.vdouble(2.1,80.00), values = cms.vdouble(0.993856055235), uncertainties = cms.vdouble(0.0071812746821,0.0071812746821)),
			cms.PSet(lowBounds = cms.vdouble(1.2,80.00), upBounds = cms.vdouble(2.1,120.00), values = cms.vdouble(0.991182209493), uncertainties = cms.vdouble(0.00753601512908,0.00753601512908)),
			cms.PSet(lowBounds = cms.vdouble(1.2,120.00), upBounds = cms.vdouble(2.1,200.00), values = cms.vdouble(0.996831513396), uncertainties = cms.vdouble(0.0090354817252,0.0090354817252)),
			cms.PSet(lowBounds = cms.vdouble(1.2,200.00), upBounds = cms.vdouble(2.1,300.00), values = cms.vdouble(0.956087938815), uncertainties = cms.vdouble(0.0183820925341,0.0183820925341)),
			cms.PSet(lowBounds = cms.vdouble(1.2,300.00), upBounds = cms.vdouble(2.1,400.00), values = cms.vdouble(1.12383072294), uncertainties = cms.vdouble(0.138534792737,0.138534792737)),
			cms.PSet(lowBounds = cms.vdouble(1.2,400.00), upBounds = cms.vdouble(2.1,800.00), values = cms.vdouble(0.903013479852), uncertainties = cms.vdouble(0.164590187976,0.164590187976)),
			cms.PSet(lowBounds = cms.vdouble(1.2,800.00), upBounds = cms.vdouble(2.1,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),

			cms.PSet(lowBounds = cms.vdouble(2.1,0.00), upBounds = cms.vdouble(2.4,20.00), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000)),
			cms.PSet(lowBounds = cms.vdouble(2.1,52.00), upBounds = cms.vdouble(2.4,55.00), values = cms.vdouble(0.916842441074), uncertainties = cms.vdouble(0.00765873327257,0.00765873327257)),			
			cms.PSet(lowBounds = cms.vdouble(2.1,55.00), upBounds = cms.vdouble(2.4,60.00), values = cms.vdouble(0.934232059345), uncertainties = cms.vdouble(0.00767503526347,0.00767503526347)),			
			cms.PSet(lowBounds = cms.vdouble(2.1,60.00), upBounds = cms.vdouble(2.4,80.00), values = cms.vdouble(0.939522854346), uncertainties = cms.vdouble(0.00762746281215,0.00762746281215)),
			cms.PSet(lowBounds = cms.vdouble(2.1,80.00), upBounds = cms.vdouble(2.4,120.00), values = cms.vdouble(0.945057510139), uncertainties = cms.vdouble(0.0100854339077,0.0100854339077)),
			cms.PSet(lowBounds = cms.vdouble(2.1,120.00), upBounds = cms.vdouble(2.4,200.00), values = cms.vdouble(0.964111458582), uncertainties = cms.vdouble(0.0221293673037,0.0221293673037)),
			cms.PSet(lowBounds = cms.vdouble(2.1,200.00), upBounds = cms.vdouble(2.4,300.00), values = cms.vdouble(0.822965708347), uncertainties = cms.vdouble(0.144571175391,0.144571175391)),
			cms.PSet(lowBounds = cms.vdouble(2.1,300.00), upBounds = cms.vdouble(2.4,400.00), values = cms.vdouble(1.00372169499), uncertainties = cms.vdouble(0.0671610296505,0.0671610296505)),
			cms.PSet(lowBounds = cms.vdouble(2.1,400.00), upBounds = cms.vdouble(2.4,800.00), values = cms.vdouble(0.749992123404), uncertainties = cms.vdouble(0.290573653738,0.290573653738)),
			cms.PSet(lowBounds = cms.vdouble(2.1,800.00), upBounds = cms.vdouble(2.4,float('inf')), values = cms.vdouble(1.0000), uncertainties = cms.vdouble(0.0000,0.0000))
		)
)


MuonIDSF = cms.PSet( MuonMethodName = cms.string("FlashggMuonWeight"),
				MethodName = cms.string("FlashggMuonFromMultiLeptonMultiJet"),
				Label = cms.string("MuonIDSF"),
				NSigmas = cms.vint32(-1,1),
				OverallRange = cms.string("abs(eta)<2.4"),
				BinList = binID,
				Debug = cms.untracked.bool(False),
				ApplyCentralValue = cms.bool(True)
				)

MuonIsoSF = cms.PSet( MuonMethodName = cms.string("FlashggMuonWeight"),
				MethodName = cms.string("FlashggMuonFromMultiLeptonMultiJet"),
				Label = cms.string("MuonIsoSF"),
				NSigmas = cms.vint32(-1,1),
				OverallRange = cms.string("abs(eta)<2.4"),
				BinList = binIso,
				Debug = cms.untracked.bool(False),
				ApplyCentralValue = cms.bool(True)
				)

MuonTriggerSF = cms.PSet( MuonMethodName = cms.string("FlashggMuonWeight"),
				MethodName = cms.string("FlashggMuonFromMultiLeptonMultiJet"),
				Label = cms.string("MuonTriggerSF"),
				NSigmas = cms.vint32(-1,1),
				OverallRange = cms.string("abs(eta)<2.4"),
				BinList = binTrigger, 
				Debug = cms.untracked.bool(False),
				ApplyCentralValue = cms.bool(True)
				)