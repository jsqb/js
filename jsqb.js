(function() {

	var jsqb = window.jsqb = function(ketString) {
		if (ketString == null && typeof ketString !== 'string') throw "[jsqb] ERROR: Invalid ket-string given. Please specify a valid ket string, e.g. '|0>', '|101>', '|0>+|1>', '|01>-|10>'.";
		var newState = new State(),
		    s = ketString.replace(/-/g,'+-').split('+');
		for (var i=0; i < s.length; i++) newState.add(new BasisState(s[i]));
		return new Qubit(newState);
	};

	// Complex class
	var Complex = jsqb.Complex = function(real, img) {
		this.real = real;
		this.img  = img || 0;
	};

	Complex.prototype.toString = function() {
		var r = Math.round(this.real * 10000) / 10000, // round to four decimal places
        i = Math.round(this.img  * 10000) / 10000;

		if (i === 0) return r+'';

		var s = i === 1 ? 'i' : (i === -1 ? '-i' : i + 'i');
		if (this.real === 0) return s;
		return r + (i < 0 ? '' : '+') + s;
	};

	Complex.add = function(v0, v1) {
		if (typeof v0 === 'number' && typeof v1 === 'number') return v0 + v1;
		if (typeof v0 === 'number') return new Complex(v1.real * v0, v1.img);
		if (typeof v1 === 'number') return new Complex(v0.real * v1, v0.img);
    return new Complex(v0.real + v1.real, v0.img + v1.img);
	};

	Complex.multiply = function(v0, v1) {
		if (typeof v0 === 'number' && typeof v1 === 'number') return v0 * v1;
		if (typeof v0 === 'number') return new Complex(v1.real * v0, v1.img * v0);
		if (typeof v1 === 'number') return new Complex(v0.real * v1, v0.img * v1);
    return new Complex(v0.real*v1.real - v0.img*v1.img, v0.real*v1.img + v0.img*v1.real);
	};

	Complex.magnitude = function(v) {
		if (typeof v === 'number') return v;
		return Math.sqrt(v.real * v.real + v.img * v.img);
	};

	// multiply matrix and vector and return a result vector
	Complex.mvmultiply = function(m, v) {
		var results = new Array(m.length);
		for(var i=0; i < m.length; i++) {
			var rowResult = 0;
			for(var j=0; j < v.length; j++) {
				rowResult = Complex.add(rowResult, Complex.multiply(m[i][j], v[j]));
			}
			results[i] = rowResult;
		}
		return results;
	};

	// STATE : a hash table to store BasisStates
	var State = jsqb.State = function(basisState) {
		this.stateMap = {}; // {stateValue: BasisState-object}
		if (basisState) this.stateMap[basisState.value] = basisState;
	};

	State.prototype = {
		add: function(newBasisState) {
			if (newBasisState == null) return this;
			var map = this.stateMap,
			    newValue = newBasisState.value,
			    oldBasisState = map[newValue];
			if (oldBasisState) {
				map[newValue] = new BasisState(Complex.add(oldBasisState.amplitude, newBasisState.amplitude), newValue, newBasisState.bitLength);
			} else {
				map[newValue] = newBasisState;
			}
			return this;
		},

		getBitLength: function() {
			var map = this.stateMap;
			for (var value in map) {
				if (map.hasOwnProperty(value)) return map[value].bitLength;
			}
			return 0; // invalid
		},

		// Return a state with fixed amplitudes following the rule, sum(amplitudes^2) = 1
		normalizeAmplitudes: function() {
			var sum = 0;
			for (var stateValue in this.stateMap) {
				var basisState = this.stateMap[stateValue],
				    magnitude = Complex.magnitude(basisState.amplitude);
						sum += magnitude * magnitude;
			}

			if (Math.abs(Complex.magnitude(sum)) > 0.0000001) {
				var scale = new Complex(1 / Math.sqrt(sum));
				for (var stateValue in this.stateMap) {
					var basisState = this.stateMap[stateValue];
					basisState.amplitude = Complex.multiply(basisState.amplitude, scale);
				}
			}
			return this;
		},

		toArray: function() {
			var basisStates = [], map = this.stateMap;
			for (var value in map) {
				if (map.hasOwnProperty(value)) basisStates.push(map[value]);
			}
			return basisStates;
		},

		clone: function() {
			var newState = new State(), map = this.stateMap;
			for (var value in map) {
				if (map.hasOwnProperty(value)) newState.add(map[value].clone());
			}
			return newState;
		},

		toString: function() {
			var buf = [], map = this.stateMap;
			for (var value in map) {
				if (map.hasOwnProperty(value)) buf.push(map[value].toString());
			}
			return 'State { ' + buf.join(', ') + ' }';
		}
	};

	// BASIS-STATE a|101> : amplitude=a, value=5, bitLength=3
	var BasisState = jsqb.BasisState = function(amplitude, value, bitLength) {
		if (typeof amplitude === 'string') { // ket string |101>, -|10>
			var ketString = amplitude.trim();
			if (ketString.indexOf('-') == 0) {
				this.amplitude = new Complex(-1);
				ketString = ketString.substring(1, ketString.length).trim();
			} else {
				this.amplitude = new Complex(1);
			}
			ketString = ketString.replace(/^\|/,'').replace(/>$/,'');
			this.value = parseInt(ketString, 2);
			this.bitLength = ketString.length;
		} else {
			if (amplitude.constructor != Complex) throw "[BasisState.constructor] ERROR: Invalid amplitude given. Please specify a Complex number.";
			this.amplitude = amplitude; // Complex number
			this.value = parseInt(value, 10);
			this.bitLength = bitLength;
		}
	};

	BasisState.prototype = {
		clone: function() {
			return new BasisState(new Complex(this.amplitude.real, this.amplitude.img), this.value, this.bitLength);
		},

		toBitString: function() {
			var binValue = this.value.toString(2);
			for (var i=0, len=this.bitLength - binValue.length; i < len; i++) {
				binValue = '0' + binValue; // zero padding
			}
			return binValue;
		},

		toString: function() {
			var a = this.amplitude.toString(); // Complex
			return a === '0' ? '' : ((a === '1' ? '' : a+" ") + "|" + this.toBitString() + ">");
		}
	};

	// MEASUREMENT : measurement result and qubit after measurement
	var Measurement = function(qubit, value) {
		this.qubit = qubit;
		this.value = value;
	};

	Measurement.prototype = {
		toString: function() {
			return this.value;
			//return "Measurement {value="+this.value+", qubit="+this.qubit+"}";
		}
	};

	// QUBIT a|01> + b|11> : state={1:a, 3:b}, bitLength=2
	var Qubit = jsqb.Qubit = function(state) {
		if (state == null) throw "[Qubit.constructor] ERROR: State object is not given. Please specify a State object.";
		this.state = state.normalizeAmplitudes();
		this.bitLength = state.getBitLength();
		if (this.bitLength == 0) throw "[Qubit.constructor] ERROR: Invalid state given. Please specify a State object containing at least one valid BasisState object.";
	};

	Qubit.prototype = {
    operate: function(srcBits, destBits, func) {
			if (func && typeof func === 'function') {
				return this._operateFunction(srcBits, destBits, func);
			}
			throw "[Qubit.operate] ERROR: Invalid function given. Please specify a function to be applied on this qubit.";
		},

    operateMatrix: function(targetBits, matrix) {
			if (matrix && matrix.length === 2 && matrix[0] && matrix[0].length === 2) {
				return this._operateGate(targetBits, function(vector) {
					return Complex.mvmultiply(matrix, vector);
				});
			}
			throw "[Qubit.operateMatrix] ERROR: Invalid matrix given. Please specify a 2x2 matrix object.";
		},

    x: function(targetBits) { // Pauli X
			return this.operateMatrix(targetBits, [[0, 1],
			                                       [1, 0]]);
		},

		hadamard: function(targetBits) {
			var SQRT1_2 = Math.SQRT1_2; // = 1.0 / Math.sqrt(2)
			return this.operateMatrix(targetBits, [[SQRT1_2,  SQRT1_2],
			                                       [SQRT1_2, -SQRT1_2]]);
		},

		tensor: function(qubit) {
			var newState = new State(),
			    thisStateMap = this.state.stateMap,
			    otherStateMap = qubit.state.stateMap;

			for (var v0 in  thisStateMap) {
				var s0 = thisStateMap[v0];
				for (var v1 in otherStateMap) {
					var s1 = otherStateMap[v1],
					    newValue = (s0.value << s1.bitLength) + s1.value,
					    newAmplitude = Complex.multiply(s0.amplitude, s1.amplitude);
					newState.add(new BasisState(newAmplitude, newValue, s0.bitLength + s1.bitLength));
				}
			}

			this.state = newState.normalizeAmplitudes();
			return this;
		},

		_operateFunction: function(srcBitRange, destBitRange, callback) {
			srcBitRange = jsqb._convertToRange(srcBitRange, this.bitLength, "[Qubit.operate] ERROR: Invalid srcBitRange (2-length array) given. Please specify 2-length array (start and end indices) or null (ALL).");
			destBitRange = jsqb._convertToRange(destBitRange, this.bitLength, "[Qubit.operate] ERROR: Invalid destBitRange (2-length array) given. Please specify 2-length array (start and end indices) or null (ALL).");

			if (srcBitRange[1] >= destBitRange[0] && destBitRange[1] >= srcBitRange[0]) {
				throw "[Qubit.operate] ERROR: srcBitRange and destBitRange cannot overlap.";
			}

			var newState = new State(),
			    srcBitMask  = ((1 << (1 + srcBitRange[1]  - srcBitRange[0]))  - 1) << srcBitRange[0];
			    destBitMask = ((1 << (1 + destBitRange[1] - destBitRange[0])) - 1) << destBitRange[0];

			for (var stateValue in this.state.stateMap) {
				var basisState = this.state.stateMap[stateValue],
				    destBitValue = stateValue & destBitMask,
				    x = (stateValue & srcBitMask) >> srcBitRange[0],
				    result = (callback(x) << destBitRange[0]) & destBitMask;
				newState.add(new BasisState(basisState.amplitude, stateValue ^ result, basisState.bitLength));
			}

			this.state = newState.normalizeAmplitudes();
			return this;
		},

		_operateGate: function(srcBits, callback) {
			srcBits = jsqb._convertToArray(srcBits, this.bitLength, "[Qubit.operate] ERROR: Invalid srcBits (index-number of bits) given. Please specify a number or an array of numbers or null (ALL).");

			for (var i=0; i < srcBits.length; i++) {
				this._operateGateWithOneBit(srcBits[i], callback);
			}
			return this;
		},

		_operateGateWithOneBit: function(srcBit, callback) {
			srcBit = srcBit || 0;

			var newState = new State(),
			    srcMask  = 1 << srcBit;

			for (var stateValue in this.state.stateMap) {
				var basisState = this.state.stateMap[stateValue],
				    srcMaskedValue  = stateValue & srcMask,
				    srcBitValue  = srcMaskedValue  != 0 ? 1 : 0;
		
				var newValues = [stateValue - srcMaskedValue, stateValue - srcMaskedValue + srcMask];
				var vector = callback([1-srcBitValue, srcBitValue], srcBitValue);

				for (var j=0; j<newValues.length; j++) {
					if (Math.abs(Complex.magnitude(vector[j])) > 0.0000001) {
						newState.add(new BasisState(Complex.multiply(basisState.amplitude, vector[j]), newValues[j], basisState.bitLength));
					}
				}
			}

			this.state = newState.normalizeAmplitudes();
		},

		clone: function() {
			return new Qubit(this.state.clone());
		},

		measure: function(srcBits) {
			srcBits = jsqb._convertToArray(srcBits, this.bitLength, "[Qubit.measure] ERROR: Invalid srcBits (index-number of bits) given. Please specify a number or an array of numbers or null (ALL).");

			var newState = this.state,
			    measurement = new Measurement(this, 0);


			for (var i=0; i < srcBits.length; i++) {
				var result = this._measureOnOneBit(newState, srcBits[i]);
				newState = result.state;
				measurement.value += result.value << i;
			}

			this.state = newState.normalizeAmplitudes();
			return measurement;
    },

		_measureOnOneBit: function(state, srcBit) {
			var newState = new State(),
			    srcMask = 1 << (srcBit || 0),
			    chosenBasisState = this._chooseBasisStateRandomly(state),
			    chosenBitValue = (chosenBasisState.value & srcMask) != 0 ? 1 : 0;

			for (var stateValue in state.stateMap) {
				var srcBitValue = (stateValue & srcMask) != 0 ? 1 : 0;
				if (srcBitValue == chosenBitValue) newState.add(state.stateMap[stateValue]);
			}

			return {state: newState, value: chosenBitValue};
    },

		_chooseBasisStateRandomly: function(state) {
			var chosenBasisState = null,
			    sum = 0,
			    threshold = Math.random(),
			    prevBasisState = null;

			for (var stateValue in state.stateMap) {
				var basisState = state.stateMap[stateValue],
				    magnitude = Complex.magnitude(basisState.amplitude);
				sum += magnitude * magnitude;
				if (sum > threshold) {
					chosenBasisState = basisState;
					break;
				}
				chosenBasisState = basisState; // set previous basis-state for magnitude errors
			}
			return chosenBasisState;
		},

		toString: function() {
			var results = [],
			    stateArray = this.state.toArray();

			stateArray.sort(function(basisState1, basisState2) {
				return basisState1.value - basisState2.value;
			});

			for (var i=0; i < stateArray.length; i++) {
				var s = stateArray[i].toString();
				if (s.length === 0) continue;

				if (results.length > 0) {
					if (s.indexOf('-') == 0) results.push(' ');
					else results.push(' +');
				}

				results.push(s);
			}
			return results.join('');
		}
	};

	// Return an array object containing arithmetic progressions.
	// start and step are optional
	jsqb.range = function(start, stop, step) {
		step = step || 1;
		if (stop == null) {
			stop = start;
			start = 0;
		}

		var nums = [];
		if (step > 0) {
			for (var i=start; i <= stop; i+=step) nums.push(i);
		} else {
			for (var i=start; i >= stop; i+=step) nums.push(i);
		}
		return nums;
	};

	jsqb._convertToArray = function(a, length, errorMessage) {
			if (a == null) {
				return jsqb.range(length-1);
			} else if (typeof a === 'number') {
				return [a];
			} else if (a.length == null) {
				throw errorMessage;
			} else {
				a.sort(function(n1, n2) {
					return n1 - n2;
				});
			}
			return a;
	};

	jsqb._convertToRange = function(a, length, errorMessage) {
			if (a == null) {
				return [0, length-1];
			} else if (typeof a === 'number') {
				return [a, a];
			} else if (a.length == null && a.length != 2) {
				throw errorMessage;
			} else {
				a.sort(function(n1, n2) {
					return n1 - n2;
				});
			}
			return a;
	};

})();

var debug = function(message, isHtml) {
	var d = document.createElement("div");
	if (isHtml === true) {
		d.innerHTML = message;
	} else {
		d.appendChild(document.createTextNode(message));
	}
	document.body.appendChild(d);
};
