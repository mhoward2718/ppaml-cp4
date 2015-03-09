function TVDScoreHMM(numStates, numSteps, groundPath, inputPath, outputPath)

    %numStates = 5;
    %numSteps = 20;
    %%  load ground marginal distribution
    %ground = csvread('problem-3-solution-smoothing.csv');
    ground = csvread(groundPath);
    %% load the predicted distribution
    %input = csvread('problem-3-solution-smoothing-test.csv');
    input = csvread(inputPath);
    %% define the output file
    %outputPath = 'problem-3-output-tvd-score.csv';
    %% verify that the ground and input files are of same dimensions
    dim = size(ground);
    groundNumSteps = dim(1);
    groundNumStates = dim(2);
    dim = size(input);
    inputNumSteps = dim(1);
    inputNumStates = dim(2);
    
    if groundNumSteps ~= inputNumSteps
        disp('Number of Steps in Ground and Input do not match! EXIT.');
        return;
    end
    if groundNumStates ~= inputNumStates
        disp('Number of States in Ground and Input do not match! EXIT.');
        return;
    end

    %% compute total variation (absolute difference) between ground and input
    diff = abs(ground - input);
    diff = sum(diff,2);
    meanDiff = mean(diff);
    varDiff = var(diff,1);
    
    %% print the scores in CSV
    diff(numSteps+1) = meanDiff;
    diff(numSteps+2) = varDiff;
    ftable = array2table(diff, 'VariableNames',{'TVD'},'RowNames',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20', 'mean', 'variance'}); 
    ftable.Properties.DimensionNames = {'TimeStep' 'TVD'};
    writetable(ftable, outputPath, 'WriteRowNames', true);
end
