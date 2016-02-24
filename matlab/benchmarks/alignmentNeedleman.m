seq1=[3:8 10:15 19:24 30:35];
seq2=[3:8 10:15  30:35] + 10;
seq1n = zeros(max(seq1),1);
for i=1:length(seq1)
    seq1n(seq1(i)) = 1;
end
seq1=seq1n;
seq2n = zeros(max(seq2),1);
for i=1:length(seq2)
    seq2n(seq2(i)) = 1;
end
seq2=seq2n;

subplot(2,1,1)
plot(seq1)
subplot(2,1,2)
plot(seq2)

% https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm

g = 0;
s = 1;

n = length(seq1);
m = length(seq2);
A = zeros(n+1,m+1);
for i=1:n
    A(i,1) = i*g;
end
for j=1:m
    A(1,j) = j*g;
end
for i=2:length(seq1)
    for j=2:length(seq2)
        match = A(i-1,j-1);
        delete = A(i-1,j);
        insert = A(i,j-1);
        if seq1(i) == seq2(j) 
            match = match + s;
        elseif seq1(i) == 1 && seq2(j) == 0
            delete = delete + g;
        elseif seq1(i) == 0 && seq2(j) == 1
            insert = insert + g;
        end
        A(i,j) = max([match,delete,insert]);
    end
end

imagesc(A)

AA = [];
AB = [];
i = length(seq1);
j = length(seq2);
while (i > 1 || j > 1) && j > 0 && i > 0
    if (i > 1 && j > 1 && A(i,j) == A(i-1,j-1) + s)
        AA = [seq1(i); AA];
        AB = [seq2(j); AB];
        i = i - 1;
        j = j - 1;
    elseif (i > 1 && j > 1 && A(i,j) == A(i-1,j-1) + g)
        AA = [seq1(i); AA];
        AB = [0; AB];
        i = i - 1;
    else
        AA = [0; AA];
        AB = [seq2(j); AB];
        j = j - 1;
    end
end

subplot(2,1,1)
plot(1:length(seq1),seq1,1:length(seq2),seq2*.5)
subplot(2,1,2)
plot(1:length(AA),AA,1:length(AB),AB*.5)
