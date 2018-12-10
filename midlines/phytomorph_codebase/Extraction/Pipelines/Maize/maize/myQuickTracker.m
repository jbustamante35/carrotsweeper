function [] = myQuickTracker(I,P)
    featureMapFunction = cornerFeatureMap(5); 
    pointExtractor = simplePointExtractor(20,.00001);
    cornerTracker = featureTrack();
    cornerTracker.setWH(30,30);
    cornerTracker.attachFeatureFunction(featureMapFunction);
    cornerTracker.attachPointSelector(pointExtractor);
    cornerTracker.track(I,P);
end