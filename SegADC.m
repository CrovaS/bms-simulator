% Settings
iteration = 5000;
bit_tot = 10;
bit_bi = 7; %or 7 or 3
bit_un = bit_tot - bit_bi;
size_LSB = 1;
stdev = 0.0125*size_LSB;

% Random transistor current variation generation
tr_sum = 2^bit_tot - 1;
tr_matrix = normrnd(size_LSB, stdev, iteration, tr_sum);

% Transistor grouping for every iteration
var_matrix = zeros(iteration, 2^bit_un+bit_bi-1);
for loop_num = 1:iteration
    % Given transistor value per iteration
    tr_row = tr_matrix(loop_num,:);
    % unary tr groupings: 2^bit_bi ... 2^bit_bi (X (2^bit_un-1))
    for un_trgroup = 1:(2^bit_un-1)
        trgroup = tr_row(1, ((un_trgroup-1)*2^bit_bi+1):(un_trgroup*2^bit_bi));
        var_matrix(loop_num, un_trgroup)=sum(trgroup, 'all');
    end
    % binary tr groupings: 2^(bit_bi-1), 2^(bit_bi-2), ... , 1 (X bit_bi)
    frontval = (2^bit_bi)*(2^bit_un-1);
    loopval = 0;
    for bi_trgroup = 1:bit_bi
        a = (frontval+loopval+1);
        b = (frontval+loopval+2^(bit_bi-bi_trgroup));
        trgroup = tr_row(1, (frontval+loopval+1):(frontval+loopval+2^(bit_bi-bi_trgroup)));
        loopval = loopval + 2^(bit_bi-bi_trgroup);
        var_matrix(loop_num, (2^bit_un-1)+bi_trgroup)=sum(trgroup, 'all');
    end
end

% Decoder for switch operation
decode_matrix = zeros(2^bit_tot, 2^bit_un-1+bit_bi);
for in = 1:2^bit_tot
    % Separate unary / binary bit
    separation_un = fix((in-1)/(2^bit_bi));
    separation_bi = mod((in-1),2^bit_bi);
    % Decode: unary part
    for dec_one = 2^bit_un-separation_un:2^bit_un-1
        decode_matrix(in, dec_one) = 1;
    end
    % Decode: binary part
    decode_matrix(in, 2^bit_un:2^bit_un-1+bit_bi) = flip(de2bi(separation_bi,bit_bi),2);
end

% Multyply matrixes
calc_matrix = var_matrix * transpose(decode_matrix);

% Get error
ideal_range = 0:2^bit_tot-1;
ideal_matrix = repmat(ideal_range, iteration, 1);
ideal_matrix = size_LSB * ideal_matrix;
error_matrix = calc_matrix - ideal_matrix;

% Calibration of INL values
cali_INL_matrix = error_matrix-(error_matrix(:,end)-error_matrix(:,1))/(2^bit_tot)*[0:1:(2^bit_tot-1)];
cali_DNL_matrix = transpose(diff(transpose(cali_INL_matrix)));

% Apply RMS, get INL
rms_INL_matrix = zeros(1,2^bit_tot);
for digital_num = 1:2^bit_tot
   cali_column = cali_INL_matrix(:,digital_num);
   rms_INL_matrix(1,digital_num) = rms(cali_column);
end

% Get DNL
rms_DNL_matrix = zeros(1,2^bit_tot-1);
for digital_num = 1:2^bit_tot-1
   cali_column = cali_DNL_matrix(:,digital_num);
   rms_DNL_matrix(1,digital_num) = rms(cali_column);
end

% Plot DNL / INL
plot(cali_INL_matrix);
figure();
plot(cali_DNL_matrix);