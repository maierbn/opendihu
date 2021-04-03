% display matrix entries for different processes, 
% fill in the mesh sizes on the processes

% data for 4 processes
%n_local_dofs_muscle = [555,555,555,555];
%n_local_dofs_fat    = [444,1036,444,1036];
n_local_dofs_muscle = [570,570,570,570];
n_local_dofs_fat    = [456,1102,459,1108];
n_columns = 5;
colors = [0.5 0 0;   1 0.8 0;     1 0.5 0.055;  0.839 0.153 0.157];
matrix = Mat_0x556060db4990_0_singleMat;

% data for 16 processes
%n_local_dofs_muscle = [190,95,95,190,190,95,95,190,180,90,90,180,180,90,90,180];
%n_local_dofs_fat    = [152,76,76,456,152,76,76,456,144,72,72,432,144,72,72,432];
%colors = [0.5,0.0,0.0;1.0,0.80000000000000004,0.0;1.0,0.49803921568627452,0.054901960784313725;			0.83921568627450982,			0.15294117647058825,			0.15686274509803921;			1.0,			1.0,			0.0;			1.0,			0.0,			1.0;			0.0,			1.0,			1.0;			0.62999923704890515,			0.62999923704890515,			1.0;			0.66999313344014644,			0.50000762951094835,			0.33000686655985351;			1.0,			0.50000762951094835,			0.74999618524452583;			0.53000686655985352,			0.3499961852445258,			0.7000076295109483;			1.0,			0.74999618524452583,			0.50000762951094835;			0.66666666666666663,0.66666666666666663,1.0;		1.0,			0.0,			0.0;			0.0,			1.0,			0.0;			0.0,			0.0,			1.0];
%matrix = Mat_0x55ebeff16510_0_singleMat;

%% computation
n_compartments = n_columns-2;
n_ranks = size(n_local_dofs_muscle,2);

% count total number of rows
n_global_dofs_muscle = sum(n_local_dofs_muscle);
n_global_dofs_fat = sum(n_local_dofs_fat);
n_rows = n_global_dofs_muscle * (n_compartments+1) + n_global_dofs_fat;
fprintf("n_rows: %d\n",n_rows)

% count number of dofs per rank
n_dofs_per_rank = zeros(1,n_ranks);
for rank_no = 1:n_ranks
   n_dofs_per_rank(rank_no) = n_local_dofs_muscle(rank_no)*(n_compartments+1) + n_local_dofs_fat(rank_no);
end

%% create original matrix and original diagonal matrix by reversing the dof reordering
original_matrix = sparse(n_rows,n_rows);
original_diagonal_matrix = sparse(n_rows,n_rows);
diagonal_matrix = sparse(n_rows,n_rows);

close all

% create first figure for the full original matrix
f1 = figure('Position',[100 100 840 800],'visible','off');
set(axes, 'Ydir', 'reverse')
set(gca,'FontSize',20)
hold on;
axis equal

% create second figure for the diagonal part of the original matrix
f2 = figure('Position',[400 100 840 800],'visible','off');
set(axes, 'Ydir', 'reverse')
set(gca,'FontSize',20)
hold on;
axis equal


% loop over row ranks
for row_rank_no = 1:n_ranks
    
    if row_rank_no == 1
        row_index_start = 1;
    else
        row_index_start = sum(n_dofs_per_rank(1:row_rank_no-1))+1;
    end
    
    % loop over compartments in that rank
    for row_compartment_no = 1:n_compartments+2

        % determine source row indices
        source_row_start = row_index_start;
        if row_compartment_no > 1
            source_row_start = row_index_start ...
                + n_local_dofs_muscle(row_rank_no)*(row_compartment_no-1);
        end
        
        if row_compartment_no < n_compartments+2
            source_row_size = n_local_dofs_muscle(row_rank_no);
        else
            source_row_size = n_local_dofs_fat(row_rank_no);
        end
        
        % source row indices are (source_row_start:source_row_start+source_row_size-1)
        
        % determine target row indices
        rank_offset = 0;
        if row_rank_no > 1
            rank_offset = sum(n_local_dofs_muscle(1:row_rank_no-1));
        end
        target_row_start = (row_compartment_no-1)*n_global_dofs_muscle + 1 + rank_offset;
        
        % loop over column ranks
        for column_rank_no = 1:n_ranks
            fprintf("row %d,%d col %d\n",row_rank_no,row_compartment_no,column_rank_no)
            
            if column_rank_no == 1
                column_index_start = 1;
            else
                column_index_start = sum(n_dofs_per_rank(1:column_rank_no-1))+1;
            end
        
            % loop over compartments in that rank
            for column_compartment_no = 1:n_compartments+2

                % determine source column indices
                source_column_start = column_index_start;
                if column_compartment_no > 1
                    source_column_start = column_index_start ...
                        + n_local_dofs_muscle(column_rank_no)*(column_compartment_no-1);
                end

                if column_compartment_no < n_compartments+2
                    source_column_size = n_local_dofs_muscle(column_rank_no);
                else
                    source_column_size = n_local_dofs_fat(column_rank_no);
                end
                
                % source column indices are (source_column_start:source_column_start+source_column_size-1)
        
                % determine target column indices
                rank_offset = 0;
                if column_rank_no > 1
                    rank_offset = sum(n_local_dofs_muscle(1:column_rank_no-1));
                end
                target_column_start = (column_compartment_no-1)*n_global_dofs_muscle + 1 + rank_offset;

                % output the current copy operation
                %fprintf("row %d,%d | column %d,%d, [%d:%d] = [%d:%d] n_rows: %d\n", ...
                %    row_rank_no, row_compartment_no, ...
                %    column_rank_no, column_compartment_no, ...
                %    target_row_start, target_row_start+source_row_size-1, ...
                %    source_row_start, source_row_start+source_row_size-1, ...
                %    n_rows)

                % read the block from the parallel reorder matrix
                source_block = matrix(...
                    source_row_start:source_row_start+source_row_size-1, ...
                    source_column_start:source_column_start+source_column_size-1 ...
                );
            
                % assign the block to the original matrix
                original_matrix(...
                    target_row_start:target_row_start+source_row_size-1, ...
                    target_column_start:target_column_start+source_column_size-1 ...
                ) = source_block;
            
                % also assign to the diagonal matrix, if we are at a
                % diagonal block
                if row_compartment_no == column_compartment_no
                    original_diagonal_matrix(...
                        target_row_start:target_row_start+source_row_size-1, ...
                        target_column_start:target_column_start+source_column_size-1 ...
                    ) = source_block;
                
                    diagonal_matrix(...
                        target_row_start:target_row_start+source_row_size-1, ...
                        target_column_start:target_column_start+source_column_size-1 ...
                    ) = source_block;
                
                    % select figure for diagonal matrix
                    %figure(f2);
                    set(0,'CurrentFigure',f2);
                
                    % create separate matrix for display
                    matrix_display = sparse(n_rows,n_rows);
                    matrix_display(...
                        target_row_start:target_row_start+source_row_size-1, ...
                        target_column_start:target_column_start+source_column_size-1 ...
                    ) = source_block;
                
                    % find nonzero entries and plot them as circles
                    [i,j,~] = find(matrix_display); 
                    color = colors(row_rank_no,:);
                    scatter(j,i,'MarkerEdgeColor',color,'Marker','.');
                end
            
                % plot sparse structure with color according to source row rank
                % select figure for full matrix
                %figure(f1);
                set(0,'CurrentFigure',f1);
                
                % create separate matrix for display
                matrix_display = sparse(n_rows,n_rows);
                matrix_display(...
                    target_row_start:target_row_start+source_row_size-1, ...
                    target_column_start:target_column_start+source_column_size-1 ...
                ) = source_block;
            
                % find nonzero entries and plot them as circles
                [i,j,~] = find(matrix_display); 
                color = colors(row_rank_no,:);
                scatter(j,i,'MarkerEdgeColor',color,'Marker','.');
            end
        end
    end
end


% add lines around the blocks
% loop over muscle blocks
for row_block_no = 1:n_compartments+2
    row_index_start = (row_block_no-1)*n_global_dofs_muscle+1;
    row_index_end = row_block_no*n_global_dofs_muscle;
        
    if row_block_no == n_compartments+2
        row_index_end = n_rows;
    end
    
    for column_block_no = 1:n_compartments+2
        column_index_start = (column_block_no-1)*n_global_dofs_muscle+1;
        column_index_end = column_block_no*n_global_dofs_muscle;
        
        if column_block_no == n_compartments+2
            column_index_end = n_rows;
        end
        
        for f = [f1 f2]
            %figure(f);
            set(0,'CurrentFigure',f);

            % horizontal lines
            line([row_index_start,row_index_end],[column_index_start,column_index_start],'Color','black');
            line([row_index_start,row_index_end],[column_index_end,column_index_end],'Color','black');

            % vertical lines
            line([row_index_start,row_index_start],[column_index_start,column_index_end],'Color','black');
            line([row_index_end,row_index_end],[column_index_start,column_index_end],'Color','black');
        end
        
    end
end

saveas(f1, 'original_matrix.png')
saveas(f2, 'original_diagonal_matrix.png')

%% show matrix by ranks 1
close all
%spy(matrix)
figure('Position',[100 100 840 800])
set(axes, 'Ydir', 'reverse')
set(gca,'FontSize',18)
%set(gca,'TickLength', [0 0])
hold on;
axis equal

n_dofs_per_rank = zeros(1,n_ranks);
for rank_no = 1:n_ranks
   n_dofs_per_rank(rank_no) = n_local_dofs_muscle(rank_no)*(n_compartments+1) + n_local_dofs_fat(rank_no);
end

% loop over matrix parts of ranks
for row_rank_no = 1:n_ranks
    if row_rank_no == 1
        row_index_start = 1;
    else
        row_index_start = sum(n_dofs_per_rank(1:row_rank_no-1))+1;
    end
    row_index_end = sum(n_dofs_per_rank(1:row_rank_no));
        
    matrix_display = matrix;

    if row_rank_no > 1
        % top
        matrix_display(1:row_index_start-1,1:end) = 0;
    end

    % bottom
    matrix_display(row_index_end+1:end,1:end) = 0;

    %spy(matrix_display);
    [i,j,s] = find(matrix_display); 
    color = colors(row_rank_no,:);
    disp(color)
    scatter(j,i,'MarkerEdgeColor',color,'Marker','.')

    % horizontal lines
    column_index_start = 1;
    column_index_end = n_rows;
    line([column_index_start,column_index_end],[row_index_start,row_index_start],'Color','black')
    line([column_index_start,column_index_end],[row_index_end,row_index_end],'Color','black')
end

saveas(gcf, 'reordered_matrix.png')

%% show diagonal matrix by ranks
close all
%spy(matrix)
figure('Position',[100 100 840 800])
set(axes, 'Ydir', 'reverse')
set(gca,'FontSize',18)
%set(gca,'TickLength', [0 0])
hold on;
axis equal

n_dofs_per_rank = zeros(1,n_ranks);
for rank_no = 1:n_ranks
   n_dofs_per_rank(rank_no) = n_local_dofs_muscle(rank_no)*(n_compartments+1) + n_local_dofs_fat(rank_no);
end

% loop over matrix parts of ranks
for row_rank_no = 1:n_ranks
    if row_rank_no == 1
        row_index_start = 1;
    else
        row_index_start = sum(n_dofs_per_rank(1:row_rank_no-1))+1;
    end
    row_index_end = sum(n_dofs_per_rank(1:row_rank_no));
        
    column_rank_no = row_rank_no;
    if column_rank_no == 1
        column_index_start = 1;
    else
        column_index_start = sum(n_dofs_per_rank(1:column_rank_no-1))+1;
    end
    column_index_end = sum(n_dofs_per_rank(1:column_rank_no));
    matrix_display = matrix;

    if row_rank_no > 1
        % left
        matrix_display(1:row_index_start-1,1:end) = 0;

        % top
        matrix_display(1:end,1:column_index_start-1) = 0;
    end

    % right
    matrix_display(row_index_end+1:end,1:end) = 0;


    % bottom
    matrix_display(1:end,column_index_end+1:end) = 0;

    % plot nonzero entries
    [i,j,s] = find(matrix_display); 
    color = colors(row_rank_no,:);
    scatter(j,i,'MarkerEdgeColor',color,'Marker','.')

    % horizontal lines
    line([row_index_start,row_index_end],[column_index_start,column_index_start],'Color','black')
    line([row_index_start,row_index_end],[column_index_end,column_index_end],'Color','black')

    % vertical lines
    line([row_index_start,row_index_start],[column_index_start,column_index_end],'Color',[0.5 0.5 0.5])
    line([row_index_end,row_index_end],[column_index_start,column_index_end],'Color',[0.5 0.5 0.5])
end

saveas(gcf, 'reordered_diagonal_matrix.png')

%% investigate eigenvalues
disp('the ten smallest eigenvalues of the matrix:')
eigs(matrix,10,'sm')

disp('the highest eigenvalues:')
eigs(matrix,10,'lm')

%% plot eigenvalues
disp('all eigenvalues:')
e1 = eigs(matrix,n_rows);
e2 = eigs(diagonal_matrix,n_rows);

figure;
plot(real(e1),'Marker','.','LineWidth',3,'Color',[0 0 0]);
box(gca,'on');
grid(gca,'on');
set(gca,'FontSize',18,'XMinorTick','on');
saveas(gcf, 'eigenvalues.png')

figure;
plot(real(e2),'Marker','.','LineWidth',3,'Color',[0 0 0]);
box(gca,'on');
grid(gca,'on');
set(gca,'FontSize',18,'XMinorTick','on');
saveas(gcf, 'eigenvalues_diagonal.png')