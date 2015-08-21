function mainloop(saveFile)
%    saveFile = 'segTest1';	
    clc;
    clearvars -except saveFile;
    close all;

    i = getenv('SLURM_NTASKS_PER_NODE');
    matlabpool(i);

    tic;                    % calculate program run time

    % helper functions to convert between frequency and wavelength
    c = 3e8;
    f2w = @(freq) c/freq;
    w2f = @(wave) c/wave;
    w2f10 = @(wave) log10(c/wave);

    dir = 'linear';
    poolobj = gcp;
    addAttachedFiles(poolobj,{[dir '/segLoc.m']});

    N =30;
    wave = 2e-6;
    freq = w2f(wave);    	% operating wavelength
    mu_r = 1;               % relative permeability
    epi_r = 1;              % relative permittivity

    % incident e field
    E_inc = 1;             
    E_vect = [0 1 0];
    % E_vect = [sin((90-120/2)*2*pi/360) cos((90-120/2)*2*pi/360) 0];  
    assump = 0;             % assume end currents are zero
    % assump = 1;             % assume end currents aren't zero

    % freqList = logspace(14.1,14.36,20); 
%    list = [6, 8,10];
    list = [5, 11, 15, 21, 37, 43, 53, 69, 79, 87, 95, 119, 125, 143, 161, 181];

    NList = zeros(size(list));
    inputImList = zeros(size(list));
    QList = zeros(size(list));
    indList = zeros(size(list));

    LList = zeros(size(list));
    aList = zeros(size(list));

    parfor i =1:length(list)
        NList(i) = list(i);    	% operating wavelength
        % analyze antenna setup
        ant1 = Antenna(dir, NList(i), freq, epi_r, mu_r, E_inc, E_vect, assump);
        inputImList(i) = ant1.inputIm;
        QList(i) = ant1.QFactor;
        indList(i) = ant1.ind;
        LList(i) = ant1.del_l*ant1.N;
        aList(i) = ant1.a;
    end

    toc

    mkdir(saveFile);
    fid = fopen([saveFile '/Summary.csv'], 'wt');
    fprintf(fid, '%s \nL, %e\na, %e\n\n', dir, LList(1), aList(1)); 
    fprintf(fid, ['# of Seg, ' 'Re(Input Im), ' 'Im(Input Im), ' 'Q, ' 'Ind (pF)\n']);
    fprintf(fid, '%d, %f, %f, %f, %f\n', [NList; real(inputImList); imag(inputImList); QList; indList*1e12]);  
    fclose(fid);
    
    qFig = figure('visible','off');
    plot(list, QList);
    % semilogx(freqList, QList);
%     xlabel('Frequency (Hz)');
    xlabel('# of Segments');
    ylabel('Q');
    print(qFig, [saveFile '/Q_seg'], '-dpng');
    

    indFig = figure('visible','off');
    plot(list, indList*1e12);
    % semilogx(freqList, indList);
%     xlabel('Frequency (Hz)');
    xlabel('# of Segments');
    ylabel('inductance (pF)');
    print(indFig, [saveFile '/ind_seg'], '-dpng');

    reInImFig = figure('visible','off');
    plot(list, real(inputImList));
    % semilogx(freqList, indList);
%     xlabel('Frequency (Hz)');
    xlabel('# of Segments');
    ylabel('Resistance (\Omega)');
    print(reInImFig, [saveFile '/Res_seg'], '-dpng');
    
    
    imInImFig = figure('visible','off');
    plot(list, imag(inputImList));
    % semilogx(freqList, indList);
%     xlabel('Frequency (Hz)');
    xlabel('# of Segments');
    ylabel('Reactance (\Omega)');
    print(imInImFig, [saveFile '/React_seg'], '-dpng');

end
