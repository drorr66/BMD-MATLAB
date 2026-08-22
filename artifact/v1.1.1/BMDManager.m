classdef BMDManager < handle
    %BMDMANAGER Canonical multiplicative BMD manager for polynomial experiments.
    %
    % This research implementation follows Rotter, Hamaguchi, Minato &
    % Yajima (1997): a polynomial x^k is encoded with Boolean-like BMD
    % variables x^(1), x^(2), x^(4), ... according to the binary digits of k.
    % Multiplication uses the paper's carry rule: (x^(2^j))^2 = x^(2^(j+1)).
    %
    % A reference is a 1x2 int64 vector [weight nodeId]. nodeId==1 is the
    % terminal function 1, so [c 1] represents the constant c.
    %
    % v0.7 intentionally fails closed if an exact integer edge-weight would
    % exceed flintmax (2^53-1). This keeps all conversions used for overflow
    % checks exact and makes canonical equality trustworthy without requiring
    % Symbolic Math Toolbox or a multiprecision dependency.

    properties (SetAccess = private)
        levels      % uint32, level 1 is x^1, level 2 is x^2, ...
        lowNode     % uint32 child node IDs
        highNode    % uint32 child node IDs
        lowWeight   % int64 child edge weights
        highWeight  % int64 child edge weights
    end

    properties (Access = private)
        uniqueTable
        addCache
        mulCache
        addHits = 0
        addMisses = 0
        mulHits = 0
        mulMisses = 0
        nodesCreated = 0
    end

    properties (Constant, Access = private)
        TERMINAL = int64(1)
        MAX_EXACT = 9007199254740991 % flintmax for double
        VAR_STRIDE = 128 % fixed bit slots per original polynomial variable
    end

    methods
        function obj = BMDManager()
            % Index 1 is the terminal node. Its arrays are placeholders.
            obj.levels = uint32(0);
            obj.lowNode = uint32(0);
            obj.highNode = uint32(0);
            obj.lowWeight = int64(0);
            obj.highWeight = int64(0);

            obj.uniqueTable = containers.Map('KeyType','char','ValueType','int64');
            obj.addCache = containers.Map('KeyType','char','ValueType','any');
            obj.mulCache = containers.Map('KeyType','char','ValueType','any');
        end

        function r = zero(obj) %#ok<MANU>
            r = int64([0 1]);
        end

        function r = one(obj) %#ok<MANU>
            r = int64([1 1]);
        end

        function r = constant(obj, c)
            c = obj.asExactInt(c);
            r = [c int64(1)];
        end

        function tf = isZero(obj, r) %#ok<MANU>
            tf = (r(1) == 0);
        end

        function tf = isConstant(obj, r) %#ok<MANU>
            tf = (r(2) == 1);
        end

        function r = scale(obj, a, c)
            c = obj.asExactInt(c);
            if c == 0 || a(1) == 0
                r = obj.zero();
                return;
            end
            r = int64([obj.checkedMul(a(1), c), a(2)]);
        end

        function r = negate(obj, a)
            r = obj.scale(a, -1);
        end

        function r = monomial(obj, exponent, coefficient)
            %MONOMIAL Return coefficient*x_1^exponent (univariate shorthand).
            if nargin < 3
                coefficient = 1;
            end
            r = obj.monomialVar(1, exponent, coefficient);
        end

        function r = monomialVar(obj, variableIndex, exponent, coefficient)
            %MONOMIALVAR Return coefficient*x_variableIndex^exponent.
            % Original variables are ordered in blocks. Each block reserves
            % VAR_STRIDE power-of-two BMD variables: x_i^1,x_i^2,x_i^4,...
            if nargin < 4
                coefficient = 1;
            end
            coefficient = obj.asExactInt(coefficient);
            if variableIndex < 1 || variableIndex ~= floor(variableIndex)
                error('BMD:BadVariable','variableIndex must be a positive integer.');
            end
            if variableIndex > floor(double(intmax('uint32'))/obj.VAR_STRIDE)
                error('BMD:BadVariable','variableIndex exceeds the v0.4 level range.');
            end
            if ~isscalar(exponent) || exponent < 0 || double(exponent) ~= floor(double(exponent))
                error('BMD:BadExponent','Exponent must be a nonnegative integer scalar.');
            end
            if double(exponent) > obj.MAX_EXACT
                error('BMD:BadExponent','v0.7 exponent must be <= flintmax.');
            end
            if coefficient == 0
                r = obj.zero();
                return;
            end
            r = obj.constant(coefficient);
            e = uint64(exponent);
            if e == 0
                return;
            end
            maxBit = floor(log2(double(e))) + 1;
            if maxBit >= obj.VAR_STRIDE
                error('BMD:ExponentBitRange','Exponent exceeds reserved variable bit block.');
            end
            for bit = maxBit:-1:1
                if bitget(e, bit) ~= 0
                    globalLevel = (variableIndex-1)*obj.VAR_STRIDE + bit;
                    r = obj.makeNode(uint32(globalLevel), obj.zero(), r);
                end
            end
        end

        function r = variable(obj, variableIndex)
            %VARIABLE Return x_variableIndex.
            r = obj.monomialVar(variableIndex,1,1);
        end

        function r = variablePow2(obj, bitLevel)
            %VARIABLEPOW2 Return x_1^(2^(bitLevel-1)).
            if bitLevel < 1 || bitLevel ~= floor(bitLevel) || bitLevel >= obj.VAR_STRIDE
                error('BMD:BadLevel','bitLevel must be in 1..VAR_STRIDE-1.');
            end
            r = obj.makeNode(uint32(bitLevel), obj.zero(), obj.one());
        end

        function r = geometricSum(obj, n)
            %GEOMETRICSUM Direct canonical builder for 1+x+...+x^n.
            %
            % v0.4 avoids the sequential-add construction used in earlier
            % experiment versions.  It recursively partitions exponents by
            % the current binary digit, so construction depth is O(log n)
            % and (for this family) every created node is reachable.
            r = obj.geometricSumShifted(n, 1);
        end

        function r = geometricSumShifted(obj, n, startBitLevel)
            %GEOMETRICSUMSHIFTED Sum_{j=0}^n x^(j*2^(startBitLevel-1)).
            % This is useful for structural-sharing experiments in which a
            % high-bit block index is independent of the low-bit payload.
            if ~isscalar(n) || n < 0 || double(n) ~= floor(double(n))
                error('BMD:BadGeometricN','n must be a nonnegative integer scalar.');
            end
            if double(n) > obj.MAX_EXACT
                error('BMD:BadGeometricN','n must be <= flintmax.');
            end
            if ~isscalar(startBitLevel) || startBitLevel < 1 || startBitLevel ~= floor(startBitLevel)
                error('BMD:BadLevel','startBitLevel must be a positive integer.');
            end
            if n > 0
                needed = floor(log2(double(n))) + 1;
                if startBitLevel + needed - 1 >= obj.VAR_STRIDE
                    error('BMD:ExponentBitRange','Shifted geometric sum exceeds reserved variable bit block.');
                end
            end

            memo = containers.Map('KeyType','char','ValueType','any');
            r = buildPart(uint64(n), uint32(startBitLevel));

            function rr = buildPart(nn, lev)
                if nn == 0
                    rr = obj.one();
                    return;
                end
                key = sprintf('%u|%.0f', lev, double(nn));
                if isKey(memo,key)
                    rr = memo(key);
                    return;
                end
                q = idivide(nn,uint64(2),'floor');
                if bitand(nn,uint64(1)) ~= 0
                    child = buildPart(q,lev+uint32(1));
                    rr = obj.makeNode(lev,child,child);
                else
                    low = buildPart(q,lev+uint32(1));
                    high = buildPart(q-1,lev+uint32(1));
                    rr = obj.makeNode(lev,low,high);
                end
                memo(key) = rr;
            end
        end

        function r = geometricGrid(obj, innerN, blocks, stridePower)
            %GEOMETRICGRID Repeated geometric blocks with shared BMD tails.
            %
            % Represents
            %   (sum_{i=0}^{innerN} x^i) *
            %   (sum_{j=0}^{blocks-1} x^(j*2^stridePower)).
            %
            % Requiring innerN < 2^stridePower makes the low-bit payload and
            % high-bit block index disjoint.  A sparse term list therefore
            % has (innerN+1)*blocks terms, while the ordered BMD can share a
            % single high-bit tail across the low-bit decision structure.
            if ~isscalar(innerN) || innerN < 0 || double(innerN) ~= floor(double(innerN))
                error('BMD:BadGrid','innerN must be a nonnegative integer.');
            end
            if ~isscalar(blocks) || blocks < 1 || double(blocks) ~= floor(double(blocks))
                error('BMD:BadGrid','blocks must be a positive integer.');
            end
            if ~isscalar(stridePower) || stridePower < 1 || stridePower ~= floor(stridePower)
                error('BMD:BadGrid','stridePower must be a positive integer.');
            end
            if stridePower >= obj.VAR_STRIDE-1
                error('BMD:BadGrid','stridePower is too large for the reserved bit block.');
            end
            stride = 2^double(stridePower);
            if double(innerN) >= stride
                error('BMD:BadGrid','Require innerN < 2^stridePower for disjoint bit ranges.');
            end
            if blocks > 1
                blockBits = floor(log2(double(blocks-1))) + 1;
                if stridePower + blockBits >= obj.VAR_STRIDE
                    error('BMD:ExponentBitRange','Grid block index exceeds reserved variable bit block.');
                end
            end

            tail = obj.geometricSumShifted(blocks-1, stridePower+1);
            memo = containers.Map('KeyType','char','ValueType','any');
            r = buildWithTail(uint64(innerN),uint32(1));

            function rr = buildWithTail(nn,lev)
                if nn == 0
                    rr = tail;
                    return;
                end
                key = sprintf('%u|%.0f',lev,double(nn));
                if isKey(memo,key)
                    rr = memo(key);
                    return;
                end
                q = idivide(nn,uint64(2),'floor');
                if bitand(nn,uint64(1)) ~= 0
                    child = buildWithTail(q,lev+uint32(1));
                    rr = obj.makeNode(lev,child,child);
                else
                    low = buildWithTail(q,lev+uint32(1));
                    high = buildWithTail(q-1,lev+uint32(1));
                    rr = obj.makeNode(lev,low,high);
                end
                memo(key) = rr;
            end
        end

        function r = indicatorExponents(obj, exponents)
            %INDICATOREXPONENTS Canonical sum of x^e for an exponent set.
            %
            % This v0.6 builder accepts an arbitrary univariate support set
            % with unit coefficients and recursively partitions exponents by
            % successive binary digits.  Identical suffix structures are
            % merged by the normal BMD unique table, so the final node count
            % directly reflects structural sharing in the support pattern.
            % No sequential polynomial additions are used.
            if isempty(exponents)
                r = obj.zero();
                return;
            end
            e = double(exponents(:).');
            if any(~isfinite(e)) || any(e < 0) || any(e ~= floor(e)) || any(e > obj.MAX_EXACT)
                error('BMD:BadExponentSet','All exponents must be exact nonnegative integers <= flintmax.');
            end
            e = unique(uint64(e),'sorted');
            if ~isempty(e)
                maxE = double(e(end));
                if maxE > 0
                    needed = floor(log2(maxE)) + 1;
                    if needed >= obj.VAR_STRIDE
                        error('BMD:ExponentBitRange','Exponent set exceeds reserved variable bit block.');
                    end
                end
            end
            r = buildSet(e,uint32(1));

            function rr = buildSet(vals,lev)
                if isempty(vals)
                    rr = obj.zero();
                    return;
                end
                if numel(vals)==1 && vals(1)==0
                    rr = obj.one();
                    return;
                end
                if double(lev) >= obj.VAR_STRIDE
                    error('BMD:ExponentBitRange','Exponent-set recursion exceeded reserved variable bit block.');
                end
                oddMask = (bitand(vals,uint64(1)) ~= 0);
                evenVals = bitshift(vals(~oddMask),-1);
                oddVals = bitshift(vals(oddMask),-1);
                low = buildSet(evenVals,lev+uint32(1));
                high = buildSet(oddVals,lev+uint32(1));
                rr = obj.makeNode(lev,low,high);
            end
        end

        function r = fromDense(obj, p)
            %FROMDENSE Build from MATLAB descending coefficient vector.
            if isempty(p)
                r = obj.zero();
                return;
            end
            p = p(:).';
            deg = numel(p)-1;
            r = obj.zero();
            for idx = 1:numel(p)
                c = obj.asExactInt(p(idx));
                if c ~= 0
                    exponent = deg-(idx-1);
                    r = obj.add(r, obj.monomial(exponent, c));
                end
            end
        end

        function r = add(obj, a, b)
            if obj.isZero(a)
                r = b;
                return;
            elseif obj.isZero(b)
                r = a;
                return;
            elseif a(2) == b(2)
                % Same canonical subgraph: addition is just edge-weight
                % addition. This is both exact and a major fast path.
                w = obj.checkedAdd(a(1), b(1));
                if w == 0
                    r = obj.zero();
                else
                    r = int64([w a(2)]);
                end
                return;
            end

            key = obj.binaryKey('A', a, b);
            if isKey(obj.addCache, key)
                obj.addHits = obj.addHits + 1;
                r = obj.addCache(key);
                return;
            end
            obj.addMisses = obj.addMisses + 1;

            top = min(obj.topLevel(a), obj.topLevel(b));
            [a0, a1] = obj.splitAt(a, top);
            [b0, b1] = obj.splitAt(b, top);
            low = obj.add(a0, b0);
            high = obj.add(a1, b1);
            r = obj.makeNode(uint32(top), low, high);
            obj.addCache(key) = r;
        end

        function r = subtract(obj, a, b)
            r = obj.add(a, obj.negate(b));
        end

        function r = multiply(obj, a, b)
            if obj.isZero(a) || obj.isZero(b)
                r = obj.zero();
                return;
            elseif obj.isConstant(a)
                r = obj.scale(b, a(1));
                return;
            elseif obj.isConstant(b)
                r = obj.scale(a, b(1));
                return;
            end

            key = obj.binaryKey('M', a, b);
            if isKey(obj.mulCache, key)
                obj.mulHits = obj.mulHits + 1;
                r = obj.mulCache(key);
                return;
            end
            obj.mulMisses = obj.mulMisses + 1;

            topA = obj.topLevel(a);
            topB = obj.topLevel(b);
            top = min(topA, topB);
            [a0, a1] = obj.splitAt(a, top);
            [b0, b1] = obj.splitAt(b, top);

            p00 = obj.multiply(a0, b0);
            p01 = obj.multiply(a0, b1);
            p10 = obj.multiply(a1, b0);

            if topA == topB
                % Paper decomposition rule 1:
                % F*G = F0G0 + TOP(F0G1+F1G0) + TOP^2 F1G1.
                p11 = obj.multiply(a1, b1);
                cross = obj.add(p01, p10);
                base = obj.makeNode(uint32(top), p00, cross);

                % TOP^2 is exactly the next power-of-two BMD variable.
                if mod(top-1,obj.VAR_STRIDE)+1 >= obj.VAR_STRIDE
                    error('BMD:ExponentBitRange','Polynomial multiplication exceeded the reserved variable bit block.');
                end
                carryVar = obj.makeNode(uint32(top + 1), obj.zero(), obj.one());
                carry = obj.multiply(carryVar, p11);
                r = obj.add(base, carry);
            else
                % Paper decomposition rule 2. One operand has no TOP term,
                % so one of a1,b1 is zero; the p11 term is kept explicitly
                % to mirror the published formula.
                p11 = obj.multiply(a1, b1);
                high = obj.add(obj.add(p01, p10), p11);
                r = obj.makeNode(uint32(top), p00, high);
            end

            obj.mulCache(key) = r;
        end

        function r = power(obj, a, n)
            if n < 0 || n ~= floor(n)
                error('BMD:BadPower','Power must be a nonnegative integer.');
            end
            result = obj.one();
            base = a;
            k = uint64(n);
            while k > 0
                if bitand(k, uint64(1)) ~= 0
                    result = obj.multiply(result, base);
                end
                k = bitshift(k, -1);
                if k > 0
                    base = obj.multiply(base, base);
                end
            end
            r = result;
        end

        function y = evaluate(obj, r, x)
            %EVALUATE Numeric evaluation. x may be scalar (x_1) or a vector
            % with one entry per original polynomial variable.
            if isempty(x)
                error('BMD:EvaluateInput','x must be a scalar or vector.');
            end
            x = double(x(:).');
            memo = containers.Map('KeyType','int64','ValueType','double');
            memo(int64(1)) = 1.0;
            y = double(r(1)) * evalNode(r(2));

            function v = evalNode(nodeId)
                if isKey(memo, nodeId)
                    v = memo(nodeId);
                    return;
                end
                idx = double(nodeId);
                lev = double(obj.levels(idx));
                [varIdx, bitIdx] = obj.decodeLevel(lev);
                if varIdx > numel(x)
                    error('BMD:EvaluateInput','Need at least %d variable values.',varIdx);
                end
                basis = x(varIdx);
                for jj = 2:bitIdx
                    basis = basis*basis;
                end
                lv = double(obj.lowWeight(idx)) * evalNode(int64(obj.lowNode(idx)));
                hv = double(obj.highWeight(idx)) * evalNode(int64(obj.highNode(idx)));
                v = lv + basis * hv;
                memo(nodeId) = v;
            end
        end

        function p = toDense(obj, r, maxDegree)
            %TODENSE Expand to MATLAB descending coefficient vector.
            % Intended for validation only, not large-scale use.
            if ~obj.isVar1Only(r)
                error('BMD:MultivariateDense','toDense is only defined for x_1-only polynomials.');
            end
            deg = obj.degree(r);
            if nargin < 3
                maxDegree = 1e6;
            end
            if deg > maxDegree
                error('BMD:DenseExpansionTooLarge', ...
                    'Degree %g exceeds maxDegree %g.', deg, maxDegree);
            end
            memo = containers.Map('KeyType','int64','ValueType','any');
            memo(int64(1)) = 1.0; % ascending coefficients for terminal 1
            asc = double(r(1)) * nodePoly(r(2));
            p = fliplr(asc);

            function coeff = nodePoly(nodeId)
                if isKey(memo, nodeId)
                    coeff = memo(nodeId);
                    return;
                end
                idx = double(nodeId);
                lev = double(obj.levels(idx));
                [~, bitIdx] = obj.decodeLevel(lev);
                shift = 2^(bitIdx-1);
                lo = double(obj.lowWeight(idx)) * nodePoly(int64(obj.lowNode(idx)));
                hi = double(obj.highWeight(idx)) * nodePoly(int64(obj.highNode(idx)));
                n = max(numel(lo), numel(hi)+shift);
                coeff = zeros(1,n);
                coeff(1:numel(lo)) = coeff(1:numel(lo)) + lo;
                hiIdx = (1:numel(hi)) + shift;
                coeff(hiIdx) = coeff(hiIdx) + hi;
                memo(nodeId) = coeff;
            end
        end

        function deg = degree(obj, r)
            if ~obj.isVar1Only(r)
                error('BMD:MultivariateDegree','degree() v0.4 is only defined for x_1-only polynomials.');
            end
            if obj.isZero(r)
                deg = -Inf;
                return;
            end
            memo = containers.Map('KeyType','int64','ValueType','double');
            memo(int64(1)) = 0.0;
            deg = nodeDegree(r(2));

            function d = nodeDegree(nodeId)
                if isKey(memo,nodeId)
                    d = memo(nodeId);
                    return;
                end
                idx = double(nodeId);
                lev = double(obj.levels(idx));
                [~, bitIdx] = obj.decodeLevel(lev);
                if obj.lowWeight(idx) == 0
                    d0 = -Inf;
                else
                    d0 = nodeDegree(int64(obj.lowNode(idx)));
                end
                if obj.highWeight(idx) == 0
                    d1 = -Inf;
                else
                    d1 = 2^(bitIdx-1) + nodeDegree(int64(obj.highNode(idx)));
                end
                d = max(d0,d1);
                memo(nodeId) = d;
            end
        end

        function n = reachableNodeCount(obj, r)
            if obj.isConstant(r)
                n = 0;
                return;
            end
            seen = false(1, numel(obj.levels));
            walk(r(2));
            n = sum(seen(2:end));

            function walk(nodeId)
                idx = double(nodeId);
                if idx == 1 || seen(idx)
                    return;
                end
                seen(idx) = true;
                walk(int64(obj.lowNode(idx)));
                walk(int64(obj.highNode(idx)));
            end
        end

        function s = stats(obj, r)
            if nargin < 2
                reachable = NaN;
                maxLev = NaN;
            else
                reachable = obj.reachableNodeCount(r);
                maxLev = obj.maxReachableLevel(r);
            end
            nInternal = numel(obj.levels)-1;
            bytesPerNode = (4+4+4+8+8);
            payloadBytes = nInternal * bytesPerNode;
            if isnan(reachable)
                reachablePayloadBytes = NaN;
                garbageNodes = NaN;
                workspaceOverReachable = NaN;
            else
                reachablePayloadBytes = reachable * bytesPerNode;
                garbageNodes = max(0,nInternal-reachable);
                workspaceOverReachable = nInternal / max(1,reachable);
            end
            s = struct( ...
                'reachable_nodes', reachable, ...
                'reachable_payload_bytes_lower_bound', reachablePayloadBytes, ...
                'total_internal_nodes', nInternal, ...
                'garbage_internal_nodes', garbageNodes, ...
                'workspace_over_reachable_ratio', workspaceOverReachable, ...
                'nodes_created', obj.nodesCreated, ...
                'unique_entries', obj.uniqueTable.Count, ...
                'add_cache_entries', obj.addCache.Count, ...
                'mul_cache_entries', obj.mulCache.Count, ...
                'add_hits', obj.addHits, ...
                'add_misses', obj.addMisses, ...
                'mul_hits', obj.mulHits, ...
                'mul_misses', obj.mulMisses, ...
                'max_level', maxLev, ...
                'node_payload_bytes_lower_bound', payloadBytes);
        end

        function tf = same(obj, a, b) %#ok<INUSL>
            %SAME O(1) canonical equality for references in the same manager.
            % Canonical node IDs are manager-local; refs from different
            % BMDManager instances must not be compared with this method.
            tf = isequal(a,b);
        end

        function [newObj,newRoots] = compact(obj, roots)
            %COMPACT Copy only nodes reachable from one or more roots.
            %
            % [newObj,newRoots] = obj.compact(root)
            % [newObj,newRoots] = obj.compact([rootA; rootB; ...])
            %
            % The returned manager starts with empty computed caches and a
            % unique table containing only the union of nodes reachable from
            % the supplied roots. This is useful both for memory accounting
            % and for cold-operation benchmarks that must exclude operand
            % construction from the timed region.
            if isempty(roots)
                newObj = BMDManager();
                newRoots = zeros(0,2,'int64');
                return;
            end
            roots = int64(roots);
            if isvector(roots) && numel(roots)==2
                roots = reshape(roots,1,2);
            end
            if size(roots,2) ~= 2
                error('BMD:BadRoots','roots must be a 1x2 ref or an Nx2 ref matrix.');
            end

            newObj = BMDManager();
            newRoots = zeros(size(roots),'int64');
            memo = containers.Map('KeyType','int64','ValueType','any');
            memo(int64(1)) = newObj.one();

            for ii = 1:size(roots,1)
                if roots(ii,1) == 0
                    newRoots(ii,:) = newObj.zero();
                else
                    base = cloneNode(roots(ii,2));
                    newRoots(ii,:) = newObj.scale(base,roots(ii,1));
                end
            end

            function nr = cloneNode(nodeId)
                if isKey(memo,nodeId)
                    nr = memo(nodeId);
                    return;
                end
                idx = double(nodeId);
                lowBase = cloneNode(int64(obj.lowNode(idx)));
                highBase = cloneNode(int64(obj.highNode(idx)));
                lowRef = newObj.scale(lowBase,obj.lowWeight(idx));
                highRef = newObj.scale(highBase,obj.highWeight(idx));
                nr = newObj.makeNode(obj.levels(idx),lowRef,highRef);
                memo(nodeId) = nr;
            end
        end

        function clearComputedCaches(obj)
            obj.addCache = containers.Map('KeyType','char','ValueType','any');
            obj.mulCache = containers.Map('KeyType','char','ValueType','any');
            obj.addHits = 0; obj.addMisses = 0;
            obj.mulHits = 0; obj.mulMisses = 0;
        end
    end

    methods (Access = private)
        function r = makeNode(obj, level, low, high)
            if obj.isZero(high)
                r = low; % BMD redundant-node reduction
                return;
            end

            % Ordering invariant: child variables must be lower in the BMD,
            % i.e. have a numerically larger level.
            if ~obj.isConstant(low) && obj.topLevel(low) <= double(level)
                error('BMD:OrderingInvariant','Low child violates variable ordering.');
            end
            if ~obj.isConstant(high) && obj.topLevel(high) <= double(level)
                error('BMD:OrderingInvariant','High child violates variable ordering.');
            end

            lw = low(1); hw = high(1);
            ln = low(2); hn = high(2);
            if lw == 0, ln = obj.TERMINAL; end
            if hw == 0, hn = obj.TERMINAL; end

            g = gcd(abs(lw), abs(hw));
            if g == 0
                r = obj.zero();
                return;
            end

            % Bryant/Chen-style gcd normalization, with an additional sign
            % convention for negative coefficients: the first nonzero
            % outgoing weight is positive. For nonnegative paper examples
            % this is exactly the published gcd normalization.
            first = lw;
            if first == 0, first = hw; end
            signFactor = int64(1);
            if first < 0, signFactor = int64(-1); end
            factor = obj.checkedMul(g, signFactor);
            nlw = idivide(lw, factor, 'fix');
            nhw = idivide(hw, factor, 'fix');

            key = sprintf('%u|%d|%d|%d|%d', uint32(level), ln, nlw, hn, nhw);
            if isKey(obj.uniqueTable, key)
                nodeId = obj.uniqueTable(key);
            else
                nodeId = int64(numel(obj.levels) + 1);
                obj.levels(end+1) = uint32(level);
                obj.lowNode(end+1) = uint32(ln);
                obj.highNode(end+1) = uint32(hn);
                obj.lowWeight(end+1) = nlw;
                obj.highWeight(end+1) = nhw;
                obj.uniqueTable(key) = nodeId;
                obj.nodesCreated = obj.nodesCreated + 1;
            end
            r = int64([factor nodeId]);
        end

        function top = topLevel(obj, r)
            if r(2) == obj.TERMINAL
                top = Inf;
            else
                top = double(obj.levels(double(r(2))));
            end
        end

        function [r0, r1] = splitAt(obj, r, top)
            if obj.isConstant(r) || obj.topLevel(r) > top
                r0 = r;
                r1 = obj.zero();
                return;
            end
            if obj.topLevel(r) ~= top
                error('BMD:SplitInvariant','splitAt called below the top level.');
            end
            idx = double(r(2));
            incoming = r(1);
            r0 = int64([obj.checkedMul(incoming, obj.lowWeight(idx)), int64(obj.lowNode(idx))]);
            r1 = int64([obj.checkedMul(incoming, obj.highWeight(idx)), int64(obj.highNode(idx))]);
            if r0(1) == 0, r0(2) = obj.TERMINAL; end
            if r1(1) == 0, r1(2) = obj.TERMINAL; end
        end

        function key = binaryKey(obj, op, a, b) %#ok<INUSL>
            % Canonicalize key order because addition/multiplication commute.
            ka = sprintf('%d:%d', a(1), a(2));
            kb = sprintf('%d:%d', b(1), b(2));
            if obj.lexLessOrEqual(ka, kb)
                key = [op '|' ka '|' kb];
            else
                key = [op '|' kb '|' ka];
            end
        end

        function tf = lexLessOrEqual(obj, a, b) %#ok<INUSL>
            % Portable lexical comparison without string class dependency.
            n = min(length(a),length(b));
            tf = true;
            for k = 1:n
                if a(k) < b(k), tf = true; return; end
                if a(k) > b(k), tf = false; return; end
            end
            tf = length(a) <= length(b);
        end

        function [varIdx, bitIdx] = decodeLevel(obj, level)
            varIdx = floor((double(level)-1)/obj.VAR_STRIDE) + 1;
            bitIdx = mod(double(level)-1,obj.VAR_STRIDE) + 1;
        end

        function tf = isVar1Only(obj, r)
            if obj.isConstant(r)
                tf = true;
                return;
            end
            seen = false(1,numel(obj.levels));
            tf = true;
            walk(r(2));
            function walk(nodeId)
                idx = double(nodeId);
                if ~tf || idx == 1 || seen(idx), return; end
                seen(idx) = true;
                [varIdx,~] = obj.decodeLevel(double(obj.levels(idx)));
                if varIdx ~= 1
                    tf = false;
                    return;
                end
                walk(int64(obj.lowNode(idx)));
                walk(int64(obj.highNode(idx)));
            end
        end

        function n = maxReachableLevel(obj, r)
            if obj.isConstant(r)
                n = 0;
                return;
            end
            seen = false(1,numel(obj.levels));
            n = 0;
            walk(r(2));
            function walk(nodeId)
                idx = double(nodeId);
                if idx == 1 || seen(idx), return; end
                seen(idx) = true;
                n = max(n, double(obj.levels(idx)));
                walk(int64(obj.lowNode(idx)));
                walk(int64(obj.highNode(idx)));
            end
        end

        function c = asExactInt(obj, v)
            if ~isscalar(v) || ~isreal(v) || ~isfinite(double(v)) || double(v) ~= fix(double(v))
                error('BMD:ExactIntegerRequired','v0.4 requires finite integer coefficients.');
            end
            if abs(double(v)) > obj.MAX_EXACT
                error('BMD:ExactRangeExceeded', ...
                    'Coefficient exceeds v0.4 exact range +/-flintmax.');
            end
            c = int64(v);
        end

        function c = checkedAdd(obj, a, b)
            d = double(a) + double(b);
            if abs(d) > obj.MAX_EXACT || ~isfinite(d)
                error('BMD:ExactRangeExceeded', ...
                    'Exact edge-weight addition exceeded +/-flintmax.');
            end
            c = int64(d);
        end

        function c = checkedMul(obj, a, b)
            d = double(a) * double(b);
            if abs(d) > obj.MAX_EXACT || ~isfinite(d)
                error('BMD:ExactRangeExceeded', ...
                    'Exact edge-weight multiplication exceeded +/-flintmax.');
            end
            c = int64(d);
        end
    end
end
