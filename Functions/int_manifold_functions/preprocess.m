function [new_data,selected_subs,group_data] = preprocess(data,labels,tsmin,tsmax,surrogate)
%Preprocess_TMH This function preprocesses the fMRI data.
%   INPUT:
%           * data: cell(1,num_subjects) containing fMRI time series
%               (   size(data{1,subject} = (num_rois,num_TRs)  )
%           * labels: cell(1,num_subjects) containing instantaneous sleep stage
%               (   size(labels{1,subject} = (num_TRs , 1)  )
%           * tsmin: scalar defining the minimum number of time points for
%               each stage to include a subject.
%           * tsmax: scalar defining the maximum number of time points
%               considered for each stage.
%           * surrogate: boolean {true: shuffle data, false: don't}.

if nargin==4
    surrogate = false;
end

% TR of the fMRI acquisition
delta_t = 2.08;                 % sampling interval

% Filtering params.
flp = .04;                      % lowpass frequency of filter
fhi = .07;                      % highpass
k   = 2;                        % 2nd order butterworth filter
fnq = 1/(2*delta_t);            % Nyquist frequency
Wn  = [flp/fnq fhi/fnq];        % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2] = butter(k,Wn); % construct the filter

% Number of brain regions
N = size(data{1,1},1);

% Var. init.
count = 0;
selected_subs = [];
timeserie_group = [];
subj_vector = [];
phases_group = [];
stage_vector = [];
pattern_vec = [];
for nsub = 1:length(data)
    % Check if there is enought time points for each sleep stage
    labels_ind = labels{nsub};
    if sum(labels_ind == 0)>= tsmin
        if sum(labels_ind == -2)>= tsmin
            if sum(labels_ind == -3)>= tsmin
                if sum(labels_ind == -4)>= tsmin
                    count = count + 1;
                    
                    new_data{1,count}.subject = nsub;
                    selected_subs = [selected_subs,nsub];
                    
                    clear timeseriedata x;
                    xs= data{1,nsub};
                    
                    if surrogate
                        % perm_inds = randperm(length(xs));
                        % xs = xs(:,perm_inds);
                        xs = phaseRandomize(xs);
                    end
                    
                    % Demean and filter the data, and get phases
                    for seed=1:N
                        x(seed,:)=demean(detrend(xs(seed,:)));
                        timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x(seed,:));    % zero phase filter the data
                        analytic_signal = hilbert(demean(timeseriedata(seed,:)));                        
                        phases(seed,:) = angle(analytic_signal);
                    end
                    
                    % Time indices of each stage
                    TimeA{count}=find(labels{1,nsub}==0);
                    Time1{count}=find(labels{1,nsub}==-2);
                    Time2{count}=find(labels{1,nsub}==-3);
                    Time3{count}=find(labels{1,nsub}==-4);
                    
                    % Select (tsmax-tsmin) TRs from the middle
                    conditions = ['A','1','2','3'];
                    for ncond = 1:4
                        time_cond = length(eval(sprintf('Time%s{count}',conditions(ncond))));
                        if time_cond > tsmax
                            aux1 = ceil((time_cond-tsmax)/2);
                            aux2 = aux1+tsmax-1;
                            eval(sprintf('Time%s{count}=Time%s{count}(aux1:aux2);',conditions(ncond),conditions(ncond)));
                        end
                    end
                    
                    % Save the indices of the selected TRs in a vector
                    inds = [ TimeA{count}; Time1{count}; Time2{count}; Time3{count}];
                    new_data{1,count}.stages = labels{1,nsub}(inds);
                    new_data{1,count}.ts = timeseriedata(:,inds);
                    new_data{1,count}.phases = phases(:,inds);
                    
                    
                    % Append group data
                    timeserie_group = horzcat(timeserie_group,timeseriedata(:,inds));    
                    subj_vector = horzcat(subj_vector,nsub*ones(1,size(timeseriedata(:,inds),2)));
                    phases_group = horzcat(phases_group,phases(:,inds));
                    stage_vector = horzcat(stage_vector,labels{1,nsub}(inds)');
                    
                     
                    % THIS PART IS NEW FOR TMH WITH ROTATION
                    % Compute the pattern now to avoid computing it in the
                    % loop. 
                    [N,T] = size(new_data{1,count}.phases);
                    Isubdiag = find(tril(ones(N),-1));
                    pattern = zeros(T,numel(Isubdiag));
  
                    for t=1:T
                        % Compute phase synchrony between each pair of brain areas at each
                        % time point
                        patt = zeros(N,N);
                        for i=1:N
                            for j=1:i-1
                                patt(i,j)=cos(adif(new_data{1,count}.phases(i,t),new_data{1,count}.phases(j,t)));
                            end
                        end
                        pattern(t,:)=patt(Isubdiag);
                    end

                    new_data{1,count}.pattern = pattern;   
                    pattern_vec = horzcat(pattern_vec,pattern');
                end
            end
        end
    end
end

% Store group data
group_data.ts           = timeserie_group;
group_data.subj_vector  = subj_vector;
group_data.phases       = phases_group;
group_data.stages       = stage_vector;
% new
group_data.pattern          = pattern_vec;
end

