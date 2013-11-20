clear all;
createConfiguration;
% generate .txt binary feature file
generate_idsfile_test(testpath);
generate_idsfile_fg_test(testpath_fg);

outputBinaryFeature(0, 'test', 'notflipped', testfilepath, testfileflippedpath, testfilepath_fg, testfileflippedpath_fg, testfilegt);
%outputBinaryFeature(0, 'test', 'flipped', testfilepath, testfileflippedpath, testfilepath_fg, testfileflippedpath_fg, testfilegt);


%train random forest
evaluation
