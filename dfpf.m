function dfpf(matrix_in, l_x, l_y)
    grid_data = reshape(str2double(regexp(matrix_in,'[+-]?\d+\.?\d*','match')),l_y,[]);
    out = grid_data;    %set grid to a new grid
    err = 1;            %error tolerance
    conv_v =1d-6;       %1.0000e-06
    k = 1;              %initialise iteration
    
    %display headings with this
    fprintf('    k |')
    for rs = 1 : l_y - 2
        for cs = 1 : l_x - 2
            fprintf('   u(%1i,%1i) |',rs,cs)
        end
    end
    fprintf('   Error   |\n')
    %display value k to be first iteration
    fprintf(' %4i |',k)
    
    %calculate the first odd sets
    out(3,3) = 0.25*(out(5,5) + out(1,1) + out(5,1) +out(1,5));
    out(2,2) = 0.25*(out(3,3) + out(1,3) + out(3,1) +out(1,1));
    out(2,4) = 0.25*(out(3,3) + out(1,3) + out(3,5) +out(1,5));
    out(4,4) = 0.25*(out(5,5) + out(3,3) + out(5,3) +out(3,5));
    out(4,2) = 0.25*(out(5,3) + out(5,1) + out(3,1) +out(3,3));
    
    %calculate the first evens sets
    out(2,3) = 0.25*(out(1,3) + out(2,2) + out(3,3) +out(2,4));
    out(3,2) = 0.25*(out(3,1) + out(2,2) + out(3,3) +out(4,2));
    out(4,3) = 0.25*(out(5,3) + out(4,4) + out(3,3) +out(4,2));
    out(3,4) = 0.25*(out(3,5) + out(2,4) + out(3,3) +out(4,4));
    
    %display calculated odd values
    for rs =l_y - 1 :-1: 2
        for cs = 2 : l_x - 1
           fprintf(' %8.6f |',out(rs,cs))
        end
    end
    err = sqrt(sum(sum((out - grid_data).^2)));
    fprintf(' %8.6f  |\n',err);
    
    grid_data = out;
    %disp(grid_data)
    
    %iterate for the result
    while err > conv_v
        k = k + 1;      %increment the loop
        fprintf(' %4i |',k)
        
        %main calculations
        for rs =l_y - 1:-1: 2
            for cs = 2 : l_x - 1
                %laplace with gaussiedel
                out(rs,cs)= 0.25*(out(rs-1,cs) + out(rs+1,cs) + out(rs,cs-1) + out(rs,cs+1));
                %laplace with jacobian
                % out(rs,cs)= 0.25*(grid_data(rs-1,cs) + grid_data(rs+1,cs) + grid_data(rs,cs-1) + grid_data(rs,cs+1));
                fprintf(' %8.6f |',out(rs,cs))
            end
        end
        %estimate error at each calculations
        err = sqrt(sum(sum((out - grid_data).^2)));
        fprintf(' %8.6f  |\n',err);
        grid_data = out;
    end
    %display the new matrix
    fprintf(' The Final Grid System \n');
    disp(grid_data);
end