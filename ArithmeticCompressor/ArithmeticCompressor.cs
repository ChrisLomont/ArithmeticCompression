/*
MIT License

Copyright (c) 2008-2021 Chris Lomont

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*	Arithmetic compression in C#
	A simple, single-file C# arithmetic compressor/decompressor

	Source version 0.6, Dec, 2021, Chris Lomont
	C# 10.0, .NET 6.0 core
*/

using System.Diagnostics;

/* Sample Program:
 
// 0. Add a using to make names smaller
using Lomont.Compression.Arithmetic;
// 1. Get some data to compress as a byte array
byte[] data = new byte[1000000]; // compress a million zeros. Could use File.ReadBytes(...)
// 2. compress the data to another byte array
List<byte> packedData = Compressor.CompressBytes(data);
// 3. decompress the data as desired
List<byte> unpackedData = Decompressor.DecompressBytes(packedData);
// 4. view the data sizes: 1000000 -> 430 -> 1000000
Console.WriteLine("Sizes {0} -> {1} -> {2}", data.Length, packedData.Count, unpackedData.Count);
// can make classes, use custom input and output sources, can use non-byte symbol tables, can use custom models... read code
*/

/* Notes:
   Compressor and decompressor classes take optionally the count of unique symbols, 256 for bytes, can customize.
   They take optional custom source and sink bit reader and writer. See the classes.
   They take custom models. Provided are 2 order-0 models: SimpleModel which is naive and slow, and the default FastModel, based on Fenwick Trees
 */


/* TODO
 * 1. Investigate outputting a byte at a time for speed
 * 2. Try other models - using array based sorted frequency tree - see solomon compression book
 * 3. Test 16 bit symbols as well as 8 bit ones for speed
 * 4. need encoder and decoder to have a FLUSH_CONTEXT symbol used when an error is about to occur, such as when underflow is about to occur or compute length possible to encode before code breaks.
 * 5. consider rewriting with range stored as [L,L+R) where L is Lower, R is Range note then that we need R>= number of symbols read at all times.
 * 6. in model, if total count > some max value, then divide all by two?! - keeps range/total > 0
 * 7. Check two ways to compute new range - speed/versus compression tradeoff
FASTER
				ulong step = range / model.Total;
				Debug.Assert(step > 0);
				high = low + step * right - 1; // -1 for open interval
				low = low + step * left;
SLOWER
				// slightly more accurate, slightly slower
				high = low + range* right/model.Total - 1; // -1 for open interval
				low = low + range* left/model.Total;
 * DONE: 8. Streaming version
 * 9. Single step version to encode parts of data at a time
 * DONE: 10. Make model external - use interface
 * DONE: 11. Investigate "Piecewise Integer Mapping for Arithmetic Coding",Stuiver & Moffat, http://www.cs.mu.oz.au/~alistair/abstracts/sm98:dcc.html
 * DONE: 12. Investigate "An Improved Data Structure for Cumulative Probability Tables", Moffat, http://www.cs.mu.oz.au/~alistair/abstracts/mof99:spe.html
 * 13. simple cases: compress file to file, byte array to byte array
 * DONE: 14. Allow more complex compression: arbitrary symbol count, backing to user place (so in mem, or to file for streaming, etc...)
 * DONE: 15. Clean comments, fix example code
 * 16. Add a order 1 (and maybe 2?) model to demonstrate
 * 17. Review "Arithmetic Coding revealed, A guided tour from theory to praxis", Bodden, Clasen, Kneis, 2007, http://www.sable.mcgill.ca/~ebodde/pubs/sable-tr-2007-5.pdf
 * */

namespace Lomont.Compression.Arithmetic;

/// <summary>
/// Class to perform arithmetic compression
/// </summary>
public class Compressor
{
	#region Simple Interface
	/// <summary>
	/// Compress a byte sequence into a compressed byte list
	/// Provides simple in memory compression
	/// </summary>
	/// <param name="data">The byte data to compress</param>
	/// <returns>The compressed byte array</returns>
	public static List<byte> CompressBytes(IEnumerable<byte> data)
	{
		var compressor = new Compressor();
		foreach (byte b in data)
			compressor.CompressSymbol(b);
		return (compressor.CompressFinish() is BitWriter m) ? m.Data : new();
	}
	#endregion

	#region Interface

	/// <summary>
	/// Construct the compressor
	/// </summary>
	/// <param name="writer">A class implementing IBitWriter to write bits to. Null for a default one.</param>
	/// <param name="symbolCount">The number of symbols possible, 256 for bytes</param>
	/// <param name="model">A class implementing IModel that models symbol probabilities. Null for a default one.</param>
	public Compressor(
		IBitWriter? writer = null,
		int symbolCount = 256,
		IModel? model = null
		)
	{
		this.model = model ?? new FastModel((ulong)symbolCount);
		this.writer = writer ?? new BitWriter();
	}

	/// <summary>
	/// compress a symbol into the stream
	/// </summary>
	/// <param name="symbol">The symbol to compress</param>
	public void CompressSymbol(ulong symbol) => EncodeSymbolHelper(symbol);

	/// <summary>
	/// call when done to flush the stream, returns the IBitWriter
	/// that received the compressed bits.
	/// </summary>
	/// <returns></returns>
	public IBitWriter CompressFinish()
	{
		FlushEncoder();
		return writer;
	}

	/// <summary>
	/// Get optimal length in bits if perfect coding used.
	/// Used to compare this implementation with perfection.
	/// </summary>
	public ulong OptimalLength => model.OptimalLength;

	#endregion

	#region Implementation

	const int BitLength = 62;  // number of bits used - todo - analyze this and make optimal?
	const ulong MaxRange = 1UL << BitLength; // highest bit used, range is [0,1] = [0,maxRange]
	const ulong HalfRange = MaxRange >> 1;     // half of the range [0,1)
	const ulong QuarterRange = HalfRange >> 1;      //  1/4 of the range [0,1)
	const ulong ThreeQuarterRange = 3 * QuarterRange;    //  3/4 of the range [0,1)

	// leave highest bits open to prevent overflow
	ulong rangeHigh = HalfRange + HalfRange - 1; // the high value of the current range [rangeLow,rangeHigh)
	ulong rangeLow;  // the low value of the current range [rangeLow,rangeHigh)

	/// <summary>
	/// track how many underflows are unaccounted for	
	/// </summary>
	long underflow;

	/// <summary>
	/// The probability model
	/// </summary>
	readonly IModel model;

	/// <summary>
	/// This member allows writing a bit at a time
	/// </summary>
	readonly IBitWriter writer;

	/// <summary>
	/// Encode a single symbol, updating internals
	/// </summary>
	/// <param name="symbol">Symbol to encode</param>
	void EncodeSymbolHelper(ulong symbol)
	{
#if DEBUG
		checked
#endif
		{
			Debug.Assert(rangeLow < rangeHigh);
			ulong range = rangeHigh - rangeLow + 1; // +1 for open interval
			model.GetRangeFromSymbol(symbol, out ulong left, out ulong right);

			ulong step = range / model.Total;
			Debug.Assert(step > 0);
			rangeHigh = rangeLow + step * right - 1; // -1 for open interval
			rangeLow += step * left;

			model.AddSymbol(symbol); // this has to be done AFTER range lookup so decoder can follow it

			// todo - analyze loops: see if needs to be 2 in 1, and see if E3 loop needs merged
			// scaling types E1, E2, E3
			while ((rangeHigh < HalfRange) || (HalfRange <= rangeLow))
			{
				if (rangeHigh < HalfRange)
				{ // E1 type scaling
					writer.Write(0);
					while (underflow > 0) { --underflow; writer.Write(1); }
					rangeHigh = (rangeHigh << 1) + 1;
					rangeLow <<= 1;
				}
				else
				{ // E2 type scaling
					writer.Write(1);
					while (underflow > 0) { --underflow; writer.Write(0); }
					rangeHigh = ((rangeHigh - HalfRange) << 1) + 1;
					rangeLow = (rangeLow - HalfRange) << 1;
				}
			}
			while ((QuarterRange <= rangeLow) && (rangeHigh < ThreeQuarterRange))
			{ // E3 type scaling
				underflow++;   // todo - if about to overflow, need a flush context symbol?!
				rangeLow = (rangeLow - QuarterRange) << 1;
				rangeHigh = ((rangeHigh - QuarterRange) << 1) + 1;
			}
			// todo - is high - low > half here? if so, assert it
			Debug.Assert(rangeHigh - rangeLow >= QuarterRange);
		}
	}

	/// <summary>
	/// Flush final data out for encoding.
	/// </summary>
	void FlushEncoder(bool writeEof = true)
	{
		if (writeEof)
			CompressSymbol(model.EndOfFileSymbol);

		// write enough bits to finalize the location of the interval
		// interval always holds at least 1/4 of the range, so cases:
		if (rangeLow < QuarterRange) // low < quarter < half <= high
		{
			writer.Write(0); // low end of the range
			for (int i = 0; i < underflow + 1; ++i) // need a 1 and then overflow bits
				writer.Write(1);
		}
		else // low < half < quarter3 <= high
		{
			writer.Write(1); // low end of range, decoder adds 0s automatically on decode
			writer.Write(0); // we'll write this just in case decoder doesn't add 0's properly..
		}
		writer.Flush();
	}

	#endregion
}

/// <summary>
/// Class to perform arithmetic decompression
/// </summary>
public class Decompressor
{
	#region Simple interface
	/// <summary>
	/// Decompress a compressed byte sequence into a decompressed byte list
	/// Provides simple in memory decompression
	/// </summary>
	/// <param name="data">The byte data to decompress</param>
	/// <returns>The decompressed byte data</returns>
	public static List<byte> DecompressBytes(IEnumerable<byte> data)
	{
		var decompressor = new Decompressor(data);
		var output = new List<byte>();
		do
		{
			var (done, symbol) = decompressor.DecompressSymbol();
			if (done) break;
			output.Add((byte)symbol);
		} while (true);
		return output;
	}
	#endregion

	#region Interface

	/// <summary>
	/// Construct the decompressor given a stream of bytes to decompress
	/// </summary>
	/// <param name="data">A stream of bytes compressed with the Compressor</param>
	/// <param name="symbolCount">The number of symbols possible, 256 for bytes</param>
	/// <param name="model">A class implementing IModel that models symbol probabilities. Null for a default one.</param>
	public Decompressor(
		IEnumerable<byte> data,
		int symbolCount = 256,
		IModel? model = null
		) : this(new BitReader(data), symbolCount, model)
	{
	}

	/// <summary>
	/// Construct the decompressor given a IBitReader that provides the stream of bits to decompress.
	/// </summary>
	/// <param name="reader">A IBitReader that provides the bits to compress</param>
	/// <param name="symbolCount">The number of symbols possible, 256 for bytes</param>
	/// <param name="model">A class implementing IModel that models symbol probabilities. Null for a default one.</param>
	public Decompressor(
		IBitReader reader,
		int symbolCount = 256,
		IModel? model = null
		)
	{
		this.model = model ?? new FastModel((ulong)symbolCount);
		this.reader = reader;

		// get initial value in [0,1) scaled range
		for (int i = 0; i < BitLength; ++i)
			currentDecodeValue = (currentDecodeValue << 1) | this.reader.Read();
	}

	/// <summary>
	/// Get a symbol. It is only valid if Done is false.
	/// </summary>
	/// <returns></returns>
	public (bool Done, ulong Symbol) DecompressSymbol()
	{
		var symbol = DecodeSymbolHelper();
		return (symbol == model.EndOfFileSymbol, symbol);
	}

	/// <summary>
	/// Get optimal length in bits if perfect coding used.
	/// Used to compare this implementation with perfection.
	/// </summary>
	public ulong OptimalLength => model.OptimalLength;

	#endregion

	#region Implementation

	const int BitLength = 62;  // number of bits used - todo - analyze this and make optimal?
	const ulong MaxRange = 1UL << BitLength; // highest bit used, range is [0,1] = [0,maxRange]
	const ulong HalfRange = MaxRange >> 1;     // half of the range [0,1)
	const ulong QuarterRange = HalfRange >> 1;      //  1/4 of the range [0,1)
	const ulong ThreeQuarterRange = 3 * QuarterRange;    //  3/4 of the range [0,1)

	// leave highest bits open to prevent overflow
	ulong rangeHigh = HalfRange + HalfRange - 1; // the high value of the current range [rangeLow,rangeHigh)
	ulong rangeLow;  // the low value of the current range [rangeLow,rangeHigh)

	/// <summary>
	/// The current value of the the decoding state
	/// This is in [rangeLow, rangeHigh]
	/// </summary>
	ulong currentDecodeValue;

	/// <summary>
	/// The probability model
	/// </summary>
	readonly IModel model;

	/// <summary>
	/// This member allows reading a bit at a time
	/// </summary>
	readonly IBitReader reader;

	/// <summary>
	/// Decode the next symbol and returns it.
	/// Once returns EOF, do not call this anymore
	/// Also updates current value in range [0,1).
	/// </summary>
	/// <returns>The decoded symbol.</returns>
	ulong DecodeSymbolHelper()
	{
#if DEBUG
		checked
#endif
		{
			Debug.Assert(rangeLow < rangeHigh);
			ulong range = rangeHigh - rangeLow + 1;

			ulong step = range / model.Total;
			Debug.Assert(step > 0);
			Debug.Assert(currentDecodeValue >= rangeLow);
			Debug.Assert(rangeHigh >= currentDecodeValue);
			ulong value = (currentDecodeValue - rangeLow) / step; // the interval location to lookup
			ulong symbol = model.GetSymbolAndRange(value, out ulong left, out ulong right);
			rangeHigh = rangeLow + step * right - 1; // -1 for open interval
			rangeLow += step * left;

			model.AddSymbol(symbol);

			// scaling types E1, E2, E3
			while ((rangeHigh < HalfRange) || (HalfRange <= rangeLow))
			{
				if (rangeHigh < HalfRange)
				{ // E1 type scaling
					rangeHigh = (rangeHigh << 1) + 1;
					rangeLow <<= 1;
					currentDecodeValue = (currentDecodeValue << 1) | reader.Read();
				}
				else
				{ // E2 type scaling
					rangeHigh = ((rangeHigh - HalfRange) << 1) + 1;
					rangeLow = (rangeLow - HalfRange) << 1;
					currentDecodeValue = ((currentDecodeValue - HalfRange) << 1) | reader.Read();
				}
			}
			while ((QuarterRange <= rangeLow) && (rangeHigh < ThreeQuarterRange))
			{ // E3 type scaling
				rangeLow = (rangeLow - QuarterRange) << 1;
				rangeHigh = ((rangeHigh - QuarterRange) << 1) + 1;
				currentDecodeValue = ((currentDecodeValue - QuarterRange) << 1) | reader.Read();
			}
			return symbol;    // todo - can do this earlier to avoid final looping?
		}
	}

	#endregion
}

#region Models
public interface IModel
{
	/// <summary>
	/// This symbol marks the end of file.
	/// </summary>
	ulong EndOfFileSymbol { get; }

	/// <summary>
	/// The total number of symbols seen
	/// </summary>
	ulong Total { get; set; }

	/// <summary>
	/// Add a new symbol to the probability table
	/// </summary>
	/// <param name="symbol">The symbol to add</param>
	void AddSymbol(ulong symbol);

	/// <summary>
	/// Given a symbol, return the range it falls into
	/// </summary>
	/// <param name="symbol">The symbol to lookup</param>
	/// <param name="low">The low end of the range</param>
	/// <param name="high">The high end of the range</param>
	void GetRangeFromSymbol(ulong symbol, out ulong low, out ulong high);

	/// <summary>
	/// Look up the symbol with given value, and return range
	/// </summary>
	/// <param name="value"></param>
	/// <param name="low">The low end of the range</param>
	/// <param name="high">The high end of the range</param>
	/// <returns>The symbol</returns>
	ulong GetSymbolAndRange(ulong value, out ulong low, out ulong high);

	/// <summary>
	/// Based on current counts, find optimal length of bits that the data could fit into
	/// </summary>
	ulong OptimalLength { get; }

} // IModel

/// <summary>
/// This model represents a simple and straightforward modeling of 
/// data probabilities.
/// Uses a simple, straightforward model that can be used to
/// benchmark other models.
/// </summary>
public class SimpleModel : IModel
{
	#region Interface

	/// <summary>
	/// End of File marker.
	/// </summary>
	public ulong EndOfFileSymbol { get; }

	/// <summary>
	/// The total number of symbols seen
	/// </summary>
	public ulong Total { get; set; }

	/// <summary>
	/// Construct a new model with the requested number of symbols
	/// </summary>
	/// <param name="size">The number of symbols needed to be stored.</param>
	public SimpleModel(ulong size)
	{
		EndOfFileSymbol = size;
		cumulativeCount = new ulong[EndOfFileSymbol + 2];
		// initialize counts to make distinct
		for (ulong i = 0; i < (ulong)(cumulativeCount.Length - 1); ++i)
			AddSymbol(i);
	}

	/// <summary>
	/// Add a new symbol to the probability table
	/// </summary>
	/// <param name="symbol">The symbol to add</param>
	public void AddSymbol(ulong symbol)
	{
		for (ulong i = symbol + 1; i < (ulong)cumulativeCount.Length; ++i)
			cumulativeCount[i]++;
		++Total;
	}

	/// <summary>
	/// Given a symbol, return the range it falls into
	/// </summary>
	/// <param name="symbol">The symbol whose range to find</param>
	/// <param name="low">The low end of the range</param>
	/// <param name="high">The high end of the range</param>
	public void GetRangeFromSymbol(ulong symbol, out ulong low, out ulong high)
	{
		ulong index = symbol;
		low = cumulativeCount[index];
		high = cumulativeCount[index + 1];
	}

	/// <summary>
	/// Look up the symbol with given value, and return range
	/// </summary>
	/// <param name="value">The value to find</param>
	/// <param name="left">The left end of the found range</param>
	/// <param name="right">The right end of the range</param>
	/// <returns>The symbol whose range includes value</returns>
	public ulong GetSymbolAndRange(ulong value, out ulong left, out ulong right)
	{
		for (ulong i = 0; i < (ulong)cumulativeCount.Length - 1; ++i)
		{
			if ((cumulativeCount[i] <= value) && (value < cumulativeCount[i + 1]))
			{
				GetRangeFromSymbol(i, out left, out right);
				return i;
			}
		}
		// if this is reached there is an unknown error elsewhere in the process
		Debug.Assert(false, "Illegal lookup overflow!");
		left = right = UInt64.MaxValue;
		return UInt64.MaxValue;
	}

	/// <summary>
	/// based on current counts, find optimal length of bits that the data could fit into
	/// </summary>
	public ulong OptimalLength
	{
		get
		{
			double sum = 0;
			double total = (Total - (ulong)cumulativeCount.Length); // number of these symbols
			for (int i = 0; i < cumulativeCount.Length - 1; ++i)
			{
				double freq = cumulativeCount[i + 1] - cumulativeCount[i] - 1;
				if (freq > 0)
				{
					double p = freq / total;
					sum += -Math.Log(p, 2) * freq;
				}
			}
			return (ulong)(sum);
		}
	}
	#endregion

	#region Implementation
	/// <summary>
	/// cumulative counts for each symbol
	/// entry j is the count of all symbol frequencies up to (but not including) symbol j
	/// </summary>
	readonly ulong[] cumulativeCount;

	#endregion

} // SimpleModel

/// <summary>
/// This model represents modeling of data probabilities 
/// using a tree based structure to provide fast O(log n) operations
/// </summary>
public class FastModel : IModel
{
	#region Implementation

	/// <summary>
	/// End of File marker.
	/// </summary>
	public ulong EndOfFileSymbol { get; }

	/// <summary>
	/// The total number of symbols seen
	/// </summary>
	public ulong Total { get; set; }

	/// <summary>
	/// Construct a new model with the requested number of symbols
	/// </summary>
	/// <param name="size"></param>
	public FastModel(ulong size)
	{
		EndOfFileSymbol = size;
		tree = new ulong[EndOfFileSymbol + 2]; // todo - this causes a slight increase over EOF+1 symbols - rethink?
											   // initialize counts to make distinct
		for (int i = 0; i < tree.Length - 1; ++i)
			AddSymbol((ulong)i);
	}

	/// <summary>
	/// Add a new symbol to the probability table
	/// </summary>
	/// <param name="symbol">The symbol to add to the table</param>
	public void AddSymbol(ulong symbol)
	{
		long k = (long)symbol;
		for (/**/ ; k < tree.Length; k |= k + 1)
			tree[k]++; // can add d here to increment element by d total instead of 1
		++Total;
	}

	/// <summary>
	/// Given a symbol, return the range it falls into
	/// </summary>
	/// <param name="symbol">The symbol to add</param>
	/// <param name="low">The low end of the range</param>
	/// <param name="high">The high end of the range</param>
	public void GetRangeFromSymbol(ulong symbol, out ulong low, out ulong high)
	{ // this uses interesting property of this tree:
	  // the parent of the higher index node of two consecutive entries will
	  // appear as an ancestor of the lower index node. This allows computing 
	  // the difference at that shared parent to get the range, then walking
	  // back on the lower index to get the lower bound.

		ulong diff = tree[symbol];
		long b = (long)(symbol);
		long target = (b & (b + 1)) - 1; // when hit this index, subtract from diff 

		ulong sum = 0;
		b = (long)(symbol - 1);
		for (; b >= 0; b = (b & (b + 1)) - 1)
		{
			if (b == target) diff -= sum;
			sum += tree[b];
		}
		if (b == target) diff -= sum;
		low = sum;
		high = sum + diff;
	}

	/// <summary>
	/// Look up the symbol with given value, and return range
	/// </summary>
	/// <param name="value">The value whose range is to be found</param>
	/// <param name="low">The low end of the range</param>
	/// <param name="high">The high end of the range</param>
	/// <returns>The symbol for this range</returns>
	public ulong GetSymbolAndRange(ulong value, out ulong low, out ulong high)
	{
		// binary search - todo - redo using bit properties of index entries
		// this takes O(log^2 n) but should take O(log n) with better searching
		ulong n = (ulong)tree.Length;
		ulong bottom = 0;
		ulong top = n;
		while (bottom < top)
		{
			ulong mid = bottom + ((top - bottom) / 2);  // Note: not (low + high) / 2 !!
			if (Query(0, (long)mid) <= value)
				bottom = mid + 1;
			else
				top = mid;
		}
		if (bottom < n)
		{ // we found a value
			GetRangeFromSymbol(bottom, out low, out high);
			Debug.Assert((low <= value) && (value < high));
			return bottom;
		}
		// if this is reached there is an unknown error elsewhere in the process
		Debug.Assert(false, "Illegal lookup overflow!");
		low = high = UInt64.MaxValue;
		return UInt64.MaxValue;
	}

	/// <summary>
	/// Find the optimal number of bits the data so far would fit into.
	/// </summary>
	public ulong OptimalLength
	{
		get
		{
			double sum = 0;
			double total = (Total - (ulong)tree.Length + 1); // number of these symbols
			for (int i = 0; i < tree.Length - 1; ++i)
			{
				double freq = Query(i, i) - 1;
				if (freq > 0)
				{
					double p = freq / total;
					sum += -Math.Log(p, 2) * freq;
				}
			}
			return (ulong)(sum);
		}
	}
	#endregion


	#region Implementation

	// Fenwick tree code based on that at http://www.algorithmist.com/index.php/Fenwick_tree
	// In this implementation, the tree is represented by an array of 64 bit unsigned integers, don't need to rescale
	// Elements are numbered by 0, 1, ..., n-1.
	// tree[i] is sum of elements with indexes i&(i+1),i&(i+2),...,i, inclusive.
	// (this is different from what is proposed in Fenwick article.)

	// todo - add scaling of tree values by half when overflow close
	// this can be done by walking tree and dividing node by half? need to prevent any freq from becoming 0

	// Returns sum of elements with indices a..b, inclusive
	ulong Query(long a, long b)
	{
		if (a == 0)
		{
			ulong sum = 0;
			for (; b >= 0; b = (b & (b + 1)) - 1)
				sum += tree[b];
			return sum;
		}

		return Query(0, b) - Query(0, a - 1);
	}

	/// <summary>
	/// Cumulative counts for each symbol stored in Fenwick tree.
	/// Entry i contains sum of elements i&(i+1), i&(i+2),...,i
	/// See above for more details on how stored.
	/// </summary>
	readonly ulong[] tree;

	#endregion
} // FastModel
#endregion

#region BitIO

public interface IBitWriter
{
	/// <summary>
	/// Output a single symbol
	/// </summary>
	/// <param name="bit"></param>
	void Write(byte bit);

	/// <summary>
	/// Fill in final byte with 0s
	/// Call before obtaining Data
	/// </summary>
	void Flush();
}

public class BitWriter : IBitWriter
{
	#region Interface
	/// <summary>
	/// Create a class that outputs bits one at a time, MSB first
	/// </summary>
	public BitWriter()
	{
		bitPos = 0;
		encodeData = new List<byte>();
		datum = 0;
	}
	/// <summary>
	/// Output a single symbol
	/// </summary>
	/// <param name="bit"></param>
	public void Write(byte bit)
	{
		datum <<= 1;
		datum = (byte)(datum | (bit & 1));
		++bitPos;
		if (bitPos == 8)
		{ // NOTE: no need to zero datum - it gets pushed along
			encodeData.Add(datum);
			bitPos = 0;
		}
	}

	/// <summary>
	/// Fill in final byte with 0s
	/// Call before obtaining Data
	/// </summary>
	public void Flush()
	{
		while (bitPos != 0)
			Write(0);
	}

	public long BitsUsed => encodeData.Count + bitPos;

	/// <summary>
	/// Obtain the internal data 
	/// </summary>
	public List<byte> Data => encodeData;
	#endregion

	#region Implementation
	int bitPos;        // bits used in current byte
	byte datum;            // current byte being created
	readonly List<byte> encodeData; // data created
	#endregion
} // BitWriter

public interface IBitReader
{
	/// <summary>
	/// Read a single bit, appending trailing zeros as needed
	/// </summary>
	/// <returns>The next bit</returns>
	byte Read();
}
public class BitReader : IBitReader
{
	#region Interface
	/// <summary>
	/// Construct a reader to parse one bit at a time, MSB first
	/// </summary>
	/// <param name="data">The data to read from</param>
	public BitReader(IEnumerable<byte> data)
	{
		bitPos = 0;
		datum = 0;
		decodeData = data.GetEnumerator();
		decodeData.MoveNext();
		datum = decodeData.Current; // first byte to output
	}

	/// <summary>
	/// Read a single bit, appending trailing zeros as needed
	/// </summary>
	/// <returns>The next bit</returns>
	public byte Read()
	{
		byte bit = (byte)(datum >> 7);
		datum <<= 1;
		++bitPos;
		if (bitPos == 8)
		{
			if (decodeData.MoveNext())
				datum = decodeData.Current; // else allow to stay at 0
			bitPos = 0;
		}
		return bit;
	}
	#endregion

	#region Implementation
	byte datum;                   // current byte of data
	int bitPos;                   // the number of bits used from datum
	readonly IEnumerator<byte> decodeData; // points to data being decoded
	#endregion
}

#endregion

// end of file