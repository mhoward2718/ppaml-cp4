function [tvd_score] = TVDScoreAgainstGround(caseposterior, samples)
    %% load ground posterior probability of each disease configuration
    caseposterior = csvread('problem-2-case1-posterior-disease-config-prob.csv');
    %% load the samples 
    samples = csvread('example-disease-samples.csv');

    dim = size(samples);
    numberOfDiseases = dim(2);

    %% binary values
    binary_values = zeros(1,numberOfDiseases);
    for i=1:numberOfDiseases
        binary_values(i) = 2^(i-1);
    end

    %% count unique samples and compute empirical mass
    [unique_samples, ia, ic] = unique(samples,'rows');
    count_unique_samples = histc(ic, 1:size(unique_samples, 1));

    empirical_mass = count_unique_samples/sum(count_unique_samples);

    %% multiply unique samples with binary values to get the index in [0 to (2^numberOfDiseases)-1]
    % empirical mass for each config

    dim = size(ia);
    total_unique = dim(1);
    unique_index_values = zeros(total_unique,numberOfDiseases);
    unique_index = zeros(total_unique,1);

    empirical_mass_config = zeros(2^numberOfDiseases,1);

    for i=1:total_unique
        unique_index_values(i,:) = unique_samples(i,:).*binary_values;
        unique_index(i) = sum(unique_index_values(i,:));
        empirical_mass_config(unique_index(i)+ 1) = empirical_mass(i);
    end

    %% absolute difference between ground posterior probability and the empirical mass
    % sum of the absolute difference is the TVD score
    abs_diff_config = abs(caseposterior - empirical_mass_config);
    tvd_score = sum(abs_diff_config);
end
