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
    Demo program

	Source version 0.6, Dec, 2021, Chris Lomont
	C# 10.0, .NET 6.0 core
*/

using C = Lomont.Compression.Arithmetic.Compressor;
using D = Lomont.Compression.Arithmetic.Decompressor;

// files compressed successfully, those with errors, those with exceptions (usually file permission stuff)
long successes = 0, errors = 0, exceptions = 0;
// total bytes compressed, resulting total size
long totalBytes = 0, totalCompressedBytes = 0; 
// elapsed milliseconds on compression, decompression
long compressionElapsedMs = 0, decompressionElapsedMs = 0;

// files to parse::
//var testPath = @"..\..\..\..\..\..\"; // big test
//var testPath = @"..\..\..\..\..\"; // medium test
//var testPath = @"..\"; // small test
var testPath = $".";

// last file with an exception in case want to explore
string lastExceptionFile = string.Empty;

// run all files
Recurse(testPath);

// last error
Console.WriteLine($"Last exception file {lastExceptionFile}");

return; // and done

// check files same
bool Check(IReadOnlyList<byte> a, byte [] b)
{
    if (a.Count != b.Length)
        return false;
    for (var i =0; i < a.Count; i++)
    {
        if (a[i] != b[i])
            return false;
    }
    return true;
}

// run a file, set stats
void TestFile(string filename)
{
    if (!File.Exists(filename))
        return;
    var originalBytes = File.ReadAllBytes(filename);
    Console.Write($"File {Path.GetFileName(filename)}");
    var originalLength = originalBytes.Length;
    totalBytes += originalLength;

    var startCompressionTime = Environment.TickCount;
    var compressedBytes = C.CompressBytes(originalBytes);
    var endCompressionTime = Environment.TickCount;
    compressionElapsedMs += endCompressionTime - startCompressionTime;
    var compressedLength = compressedBytes.Count;
    totalCompressedBytes += compressedLength;
    var compressedKbs = compressionElapsedMs == 0 ? 0 : (1000L * totalBytes) / ((1L << 10) * compressionElapsedMs);

    var startDecompressionTime = Environment.TickCount;
    var decompressedBytes = compressedBytes.Count != 0 ? D.DecompressBytes(compressedBytes) : new();
    var endDecompressionTime = Environment.TickCount;
    decompressionElapsedMs += endDecompressionTime - startDecompressionTime;
    var decompressedKbs = decompressionElapsedMs == 0 ? 0 : (1000L * totalBytes) / ((1L << 10) * decompressionElapsedMs);

    var matches = Check(decompressedBytes,originalBytes);
    if (matches) successes++; else errors++;

    Console.WriteLine($" ({Pct(compressedLength,originalLength)}) match {matches}, success {successes}, errors {errors}, except {exceptions}, com,dec kB/s {compressedKbs}:{decompressedKbs}, total ratio {Pct(totalCompressedBytes,totalBytes)}");

    string Pct(long num, long den) => den == 0 ? "0%" : $"{100 * num / den}%";
}

// recurse on path, check files
void Recurse(string path)
{
    foreach (var f in Directory.GetFiles(path))
    {
        try
        {
            TestFile(f);
        }
        catch (Exception)
        {
            Console.WriteLine("EXCEPTION....");
            lastExceptionFile = f;
            exceptions++;
        }
    }
    foreach (var dir in Directory.GetDirectories(path))
        Recurse(dir);
}




