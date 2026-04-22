clc,clear;
% delete(gcp('nocreate'));
% parpool(6);
% values = [20,21,22,23,24,25];  % 定義 w 的值
% values = [15,16,17,18,19,20];  % 定義 w 的值
% values = [4,17,18,19,20,21];
% T_Values = [1,5,7,35]; ,3,5,11,15,33,55,165

values = [4];
T_Values = [1,12];

tic
for idx = 1:length(values)

    w = values(idx);  % 根據迴圈索引取得對應的值

    for t_idx = 1:length(T_Values)
        T = T_Values(t_idx);  % 取得當前 T 的值

        p = 5;
        r = 3;
        % 1. 讓使用者輸入想選幾個序列
        se = 4;
        q = ceil(2*se/r);
        sequences = [3,4,5,6];

        % fprintf("T = %d / w = %d", T, w);
        % ✅ **每個 worker 內部都重新定義 `output_dir`**
        output_dir = sprintf('Case2_10w_M%dR%d', se, r);  % 確保 worker 知道 output 目錄
        T_folder = fullfile(output_dir, sprintf('T_%d', T));  % T 值的資料夾

        % ✅ **在 worker 內部確認目錄是否存在**
        if exist(output_dir, 'dir') == 0
            mkdir(output_dir);
        end
        if exist(T_folder, 'dir') == 0
            mkdir(T_folder);
        end
        %tic
        %% 做出Sgset
        SgSet = cell(p+1, 1);

        for g = 0:p-1
            gSet = zeros(w, 2);
            Sg = zeros(1, p*q);

            for j = 0:w-1
                gSet(j+1, :) = [mod(g*j, p), mod(j, q)];
            end

            for t = 0:p*q-1
                St = [mod(t, p), mod(t, q)];
                if any(ismember(gSet, St, 'rows'))
                    Sg(t+1) = 1;
                else
                    Sg(t+1) = 0;
                end
            end

            SgSet{g+1} = Sg;
        end

        Stmp = zeros(1, p*q);
        gSet = zeros(w, 2);
        Sg = zeros(1, p*q);

        for j = 0:w-1
            gSet(j+1, :) = [mod(j, p), mod(0, q)];
        end

        for t = 0:p*q-1
            St = [mod(t, p), mod(t, q)];
            if any(ismember(gSet, St, 'rows'))
                Stmp(t+1) = 1;
            else
                Stmp(t+1) = 0;
            end
        end
        SgSet{p+1} = Stmp;

        %% 跑loop
        ns = 100000;
        
        ta = zeros(1, se+1);
       

        for m = 1:ns
            randomSgSet = cell(se, 1);
            indexCell = cell(se, 1);

            for g = 1:p+1
                if ismember(g, sequences)
                    random_shift = randi([0, p*q]);
                    tmp = circshift(SgSet{g}, random_shift);
                    randomSgSet{g} = tmp;
                    nbo = find(randomSgSet{g} == 1);
                    indexCell{g} = nbo-1;
                else
                    randomSgSet{g} = zeros(1, p*q);
                end
            end

            addSg = zeros(1, length(SgSet{1}));

            for i = 1:length(randomSgSet)
                addSg = addSg + randomSgSet{i};
            end

            rboth = find(addSg > 0 & addSg <= r);

            sgSuccessCounts = cell(se,1);
            groupDelayTmp = 0;

            for i = 1:se
                sgSuccessCounts{i} = intersect(rboth-1, indexCell{sequences(i)});

            end

            for i = 1:length(sgSuccessCounts)
                if ~isempty(sgSuccessCounts{i})
                    first_value = sgSuccessCounts{i}(1);


                    one_Index = getOneIndextmp(sgSuccessCounts{i}, T, p, q);
                    second_group = one_Index + p*q;
                    extended_Index = [one_Index, second_group];

                    extended_AoI = zeros(1,2*p*q);

                    for d = 1:2*p*q
                        if ismember(d,extended_Index)
                            extended_AoI(d) = mod(d,T)+ 1;
                        elseif d == 1
                            extended_AoI(d) = 0;
                        else
                            extended_AoI(d) = extended_AoI(d-1) + 1;
                        end
                    end

                    calculate_AoI = extended_AoI(p*q+1:2*p*q);
                    AoI_average = sum(calculate_AoI)/(p*q);
                    ta(i) = ta(i) + AoI_average;
                end
            end
            


        end

        tmpm = p*q*ns;

        % 打開文件
        % 儲存檔案到對應的資料夾
        filename = fullfile(T_folder, sprintf('%d結果%d.xlsx', T, w));

        
        taavg = ta/ns;
        
        taavg(length(taavg))= sum(ta) / (se*ns);
        
        mat = taavg;
        writematrix(mat, filename);
    end
end
toc

