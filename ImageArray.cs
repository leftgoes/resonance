using System.Runtime.Versioning;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;
using System.IO;


namespace Resonance
{
    [SupportedOSPlatform("windows")]
    public class ImageArray
    {
        
        public readonly int Width, Height;
        private double[,] Img;

        public ImageArray(int width, int height) 
        {
            Width = width;
            Height = height;
            Img = new double[Height, Width];
        }

        public IEnumerable<(int y, int x)> Indices()
        {
            for (int i = 0; i < Height; i++) {
                for (int j = 0; j < Width; j++) {
                    yield return (i, j);
                }
            }
        }

        public double this[int y, int x]
        {
            get => Img[y, x];
            set => Img[y, x] = value;
        }

        public int Stride { get => Width % 4 == 0 ? Width : Width + 4 - Width % 4; }

        public byte[] ToByteArray()
        {
            double arrayMax = Img.Cast<double>().Max();
            byte[] bytes = new byte[Stride * Height];
            foreach ((int y, int x) in Indices())
                bytes[y * Stride + x] = (byte)Math.Round(255 * Img[y, x] / arrayMax);

            return bytes;
        }

        public Bitmap ToBitmap() // https://stackoverflow.com/a/19449233
        {
            PixelFormat formatOutput = PixelFormat.Format8bppIndexed;
            Rectangle rect = new(0, 0, Width, Height);

            Bitmap bmp = new(Stride, Height, formatOutput);
            BitmapData bmpData = bmp.LockBits(rect, ImageLockMode.ReadOnly, formatOutput);

            byte[] bytes = ToByteArray();

            Marshal.Copy(bytes, 0, bmpData.Scan0, bytes.Length);
            bmp.UnlockBits(bmpData);

            ColorPalette palette = bmp.Palette;
            Color[] entries = palette.Entries;
            for (int i = 0; i < 256; i++)
            {
                Color b = new();
                b = Color.FromArgb((byte)i, (byte)i, (byte)i);
                entries[i] = b;
            }
            bmp.Palette = palette;

            return bmp;
        }

        public void SaveToPng(string filename)
        {
            ToBitmap().Save(filename);
        }
    }

    [SupportedOSPlatform("windows")]
    public class Videoframes
    {
        public readonly int Width;
        public readonly int Height;
        ImageArray[] Frames;

        public Videoframes(int width, int height, int framesCount)
        {
            Width = width;
            Height = height;
            Frames = new ImageArray[framesCount];
            for (int i = 0; i < framesCount; i++)
                Frames[i] = new ImageArray(width, height);
        }

        public ImageArray this[int i]
        {
            get => Frames[i];
        }

        public int FramesCount { get => Frames.Length; }

        public void DrawPoint(int frameIndex, double i, double j)
        {
            int iInt = (int)i;
            int jInt = (int)j;
            double iFrac = i - (double)iInt;
            double jFrac = j - (double)jInt;

            DrawPointInt(frameIndex, iInt, jInt, (1 - iFrac) * (1 - jFrac));
            DrawPointInt(frameIndex, iInt, jInt + 1, (1 - iFrac) * jFrac);
            DrawPointInt(frameIndex, iInt + 1, jInt, iFrac * (1 - jFrac));
            DrawPointInt(frameIndex, iInt + 1, jInt + 1, iFrac * jFrac);
        }

        public void DrawPointInt(int frameIndex, int i, int j, double val)
        {
            if (0 <= i && i < Width && 0 <= j && j < Height)
            {
                Frames[frameIndex][j, i] += val;
                Console.WriteLine("Drawing");
            }
                
        }

        public void SaveToPngs(string directory)
        {
            if (Directory.Exists(directory))
                Directory.CreateDirectory(directory);

            for (int i = 0; i < FramesCount; i++)
            {
                Console.Write($"\rWriting frame {i}/{FramesCount}");
                Frames[i].SaveToPng(Path.Join(directory, "frm" + i.ToString().PadLeft(5, '0') + ".png"));
            }
        }
        
}
}
