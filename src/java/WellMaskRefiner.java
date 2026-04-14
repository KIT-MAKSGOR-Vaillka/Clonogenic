import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import javax.imageio.ImageIO;

public final class WellMaskRefiner {
    private static final double DEFAULT_AREA_SCALE = 1.02;
    private static final double LEFT_COLUMN_AREA_SCALE = 1.005;
    private static final int SAMPLE_COUNT = 240;
    private static final int SEARCH_OFFSET = 14;
    private static final int SEARCH_STEP = 2;
    private static final double CENTER_PENALTY = 0.075;

    private record Seed(int wellIndex, double x, double y, double radius) {}

    private record Result(int wellIndex, double centerX, double centerY, double outerRadius, double maskRadius) {}

    public static void main(String[] args) throws Exception {
        if (args.length != 3) {
            throw new IllegalArgumentException("Usage: java WellMaskRefiner <image> <seed_csv> <output_csv>");
        }

        BufferedImage image = ImageIO.read(new File(args[0]));
        if (image == null) {
            throw new IOException("Could not read image: " + args[0]);
        }

        double[][] gray = toGray(image);
        double[][] smoothed = gaussian5x5(gray);
        double[][] edge = sobelMagnitude(smoothed);

        List<Seed> seeds = readSeeds(Path.of(args[1]));
        List<Result> results = new ArrayList<>();
        for (Seed seed : seeds) {
            results.add(refineSeed(edge, seed));
        }
        writeResults(Path.of(args[2]), results);
    }

    private static List<Seed> readSeeds(Path path) throws IOException {
        List<Seed> seeds = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            String line = reader.readLine();
            if (line == null) {
                return seeds;
            }
            while ((line = reader.readLine()) != null) {
                if (line.isBlank()) {
                    continue;
                }
                String[] parts = line.split(",");
                seeds.add(
                    new Seed(
                        Integer.parseInt(parts[0].trim()),
                        Double.parseDouble(parts[1].trim()),
                        Double.parseDouble(parts[2].trim()),
                        Double.parseDouble(parts[3].trim())
                    )
                );
            }
        }
        return seeds;
    }

    private static void writeResults(Path path, List<Result> results) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            writer.write("well_index,center_x,center_y,outer_radius,mask_radius\n");
            for (Result result : results) {
                writer.write(
                    String.format(
                        Locale.US,
                        "%d,%.4f,%.4f,%.4f,%.4f%n",
                        result.wellIndex,
                        result.centerX,
                        result.centerY,
                        result.outerRadius,
                        result.maskRadius
                    )
                );
            }
        }
    }

    private static Result refineSeed(double[][] edge, Seed seed) {
        double bestScore = Double.NEGATIVE_INFINITY;
        double bestX = seed.x;
        double bestY = seed.y;
        double bestRadius = seed.radius;

        int minRadius = Math.max(10, (int) Math.round(seed.radius * 0.92));
        int maxRadius = Math.max(minRadius + 2, (int) Math.round(seed.radius * 1.08));

        for (int dy = -SEARCH_OFFSET; dy <= SEARCH_OFFSET; dy += SEARCH_STEP) {
            for (int dx = -SEARCH_OFFSET; dx <= SEARCH_OFFSET; dx += SEARCH_STEP) {
                double centerX = seed.x + dx;
                double centerY = seed.y + dy;
                for (int radius = minRadius; radius <= maxRadius; radius += 2) {
                    double offsetPenalty = CENTER_PENALTY * (Math.hypot(dx, dy) / seed.radius);
                    double score = circleScore(edge, centerX, centerY, radius) + 0.003 * (radius / seed.radius) - offsetPenalty;
                    if (score > bestScore) {
                        bestScore = score;
                        bestX = centerX;
                        bestY = centerY;
                        bestRadius = radius;
                    }
                }
            }
        }

        double areaScale = ((seed.wellIndex - 1) % 3 == 0) ? LEFT_COLUMN_AREA_SCALE : DEFAULT_AREA_SCALE;
        return new Result(seed.wellIndex, bestX, bestY, bestRadius, bestRadius * areaScale);
    }

    private static double circleScore(double[][] edge, double centerX, double centerY, double radius) {
        int height = edge.length;
        int width = edge[0].length;
        double sum = 0.0;
        int count = 0;

        for (int i = 0; i < SAMPLE_COUNT; i++) {
            double angle = (2.0 * Math.PI * i) / SAMPLE_COUNT;
            double x = centerX + radius * Math.cos(angle);
            double y = centerY + radius * Math.sin(angle);
            if (x < 1.0 || y < 1.0 || x >= width - 2.0 || y >= height - 2.0) {
                continue;
            }
            sum += bilinear(edge, x, y);
            count++;
        }

        if (count == 0) {
            return Double.NEGATIVE_INFINITY;
        }
        return sum / count;
    }

    private static double bilinear(double[][] image, double x, double y) {
        int x0 = (int) Math.floor(x);
        int y0 = (int) Math.floor(y);
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        double dx = x - x0;
        double dy = y - y0;

        double v00 = image[y0][x0];
        double v10 = image[y0][x1];
        double v01 = image[y1][x0];
        double v11 = image[y1][x1];

        double top = v00 * (1.0 - dx) + v10 * dx;
        double bottom = v01 * (1.0 - dx) + v11 * dx;
        return top * (1.0 - dy) + bottom * dy;
    }

    private static double[][] toGray(BufferedImage image) {
        int width = image.getWidth();
        int height = image.getHeight();
        double[][] gray = new double[height][width];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int rgb = image.getRGB(x, y);
                int red = (rgb >> 16) & 0xFF;
                int green = (rgb >> 8) & 0xFF;
                int blue = rgb & 0xFF;
                gray[y][x] = (0.2126 * red + 0.7152 * green + 0.0722 * blue) / 255.0;
            }
        }
        return gray;
    }

    private static double[][] gaussian5x5(double[][] image) {
        double[] kernel = new double[] {1.0, 4.0, 6.0, 4.0, 1.0};
        double norm = 16.0;
        return convolveVertical(convolveHorizontal(image, kernel, norm), kernel, norm);
    }

    private static double[][] convolveHorizontal(double[][] image, double[] kernel, double norm) {
        int radius = kernel.length / 2;
        int height = image.length;
        int width = image[0].length;
        double[][] out = new double[height][width];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double sum = 0.0;
                for (int k = -radius; k <= radius; k++) {
                    int xx = clamp(x + k, 0, width - 1);
                    sum += image[y][xx] * kernel[k + radius];
                }
                out[y][x] = sum / norm;
            }
        }
        return out;
    }

    private static double[][] convolveVertical(double[][] image, double[] kernel, double norm) {
        int radius = kernel.length / 2;
        int height = image.length;
        int width = image[0].length;
        double[][] out = new double[height][width];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double sum = 0.0;
                for (int k = -radius; k <= radius; k++) {
                    int yy = clamp(y + k, 0, height - 1);
                    sum += image[yy][x] * kernel[k + radius];
                }
                out[y][x] = sum / norm;
            }
        }
        return out;
    }

    private static double[][] sobelMagnitude(double[][] image) {
        int height = image.length;
        int width = image[0].length;
        double[][] out = new double[height][width];
        for (int y = 1; y < height - 1; y++) {
            for (int x = 1; x < width - 1; x++) {
                double gx =
                    -image[y - 1][x - 1] + image[y - 1][x + 1]
                    - 2.0 * image[y][x - 1] + 2.0 * image[y][x + 1]
                    - image[y + 1][x - 1] + image[y + 1][x + 1];
                double gy =
                    image[y - 1][x - 1] + 2.0 * image[y - 1][x] + image[y - 1][x + 1]
                    - image[y + 1][x - 1] - 2.0 * image[y + 1][x] - image[y + 1][x + 1];
                out[y][x] = Math.hypot(gx, gy);
            }
        }
        return out;
    }

    private static int clamp(int value, int min, int max) {
        return Math.max(min, Math.min(max, value));
    }
}
