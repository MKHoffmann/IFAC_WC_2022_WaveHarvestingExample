function s = SimplexNoise(m,args)
    arguments
        m               (1,1) {mustBeNumeric}
        args.Seed       (1,1) {mustBeInteger} = 2
    end


    rng(args.Seed);

    s = zeros([m,m]);     % Prepare output image (size: m x m)
    w = m;
    i = 0;
    while w > 3
        i = i + 1;
        % least number of random values necessary to have size(d) >= m
        n = ceil((m-1)/2^(i-1)+1); 
        d = interp2(randn([n,n]), i-1, 'spline');
        s = s + i * d(1:m, 1:m);
        w = w - ceil(w/2 - 1);
    end
    s = (s - min(s, [], 'all')) ./ (max(s, [], 'all') - min(s, [], 'all'));
%       s = histeq(s);
end