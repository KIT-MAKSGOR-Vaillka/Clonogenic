import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.gui.Wand;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.EDM;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.frame.RoiManager;
import ij.process.AutoThresholder;
import ij.process.ByteProcessor;
import ij.process.FloatPolygon;
import ij.process.FloatProcessor;
import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import javax.imageio.ImageIO;

public final class ClonogenicAnalyzer {
    static {
        System.setProperty("java.awt.headless", "true");
    }

    private static final String IJ_JAR = "/home/maksegr/Applications/Fiji.app/jars/ij-1.54p.jar";
    private static final Locale LOCALE = Locale.US;
    private static final Font LABEL_FONT = new Font("SansSerif", Font.PLAIN, 14);
    private static final Font WELL_FONT = new Font("SansSerif", Font.BOLD, 18);

    // Параметры ручной подстройки детекции колоний на TIFF-сканах и фотографиях.
    // Насколько сильно размывается локальный фон перед поиском темных и фиолетовых колоний.
    // Увеличение: фон сглаживается сильнее, шум подавляется лучше, но слабые колонии могут потеряться.
    // Уменьшение: чувствительность к мелким деталям растет, но может полезть шум. Разумно: 0.03-0.10.
    private static final double TUNE_BACKGROUND_SIGMA_SCALE = 0.07;
    // Насколько жестко используется верхний квантиль при переводе score map в бинарную маску.
    // Увеличение: маска становится строже, мелкие и слабые кандидаты чаще пропускаются.
    // Уменьшение: маска становится мягче, но может добавляться мусор. Разумно: 0.90-0.99.
    private static final double TUNE_SCORE_QUANTILE = 0.965;

    // КЛЮЧЕВОЙ ПАРАМЕТР: МИНИМАЛЬНЫЙ АБСОЛЮТНЫЙ РАЗМЕР КОЛОНИИ В ПИКСЕЛЯХ.
    // Увеличение: мелкие точки и мусор отсекаются сильнее.
    // Уменьшение: начинает считаться больше мелких колоний и мелкого шума. Разумно: 10-30.
    private static final int TUNE_MIN_COLONY_AREA = 14;

    // Смешивает строгий и мягкий пороги.
    // Увеличение: порог строже, слабые колонии чаще теряются.
    // Уменьшение: детекция мягче, но может появляться лишний фон. Разумно: 0.60-0.90.
    private static final double TUNE_THRESHOLD_RELAXATION = 0.69;
    // Насколько большим должен быть объект, чтобы мы начали считать его возможной слипшейся колонией.
    // Увеличение: разбиение слипшихся колоний включается реже.
    // Уменьшение: разбиение включается чаще, но можно получить лишнее дробление. Разумно: 2.0-5.0.
    private static final double TUNE_CLUMP_AREA_FACTOR = 3.4;
    // Радиус поиска соседних колоний при оценке локального типичного размера колонии.
    // Увеличение: локальная оценка размера учитывает более дальних соседей.
    // Уменьшение: оценка становится более локальной и нестабильной. Разумно: 0.35-0.80.
    private static final double TUNE_NEIGHBOR_RADIUS_SCALE = 0.55;
    // На сколько пикселей сужается маска лунки внутрь, чтобы построить более безопасную внутреннюю зону подсчета.
    // Увеличение: меньше влияния бортика, но краевые колонии теряются чаще.
    // Уменьшение: ближе считаем к краю, но может полезть мусор с бортика. Разумно: 0-12.
    private static final int TUNE_COUNT_MASK_EROSION_PX = 5;
    // Ширина внешней полосы у края лунки, где артефакты отбрасываются наиболее жестко.
    // Увеличение: больше защита от артефактов края.
    // Уменьшение: больше шанс посчитать колонии у борта, но и больше риск мусора. Разумно: 2-10 px.
    private static final double TUNE_OUTER_REJECT_BAND_PX = 4.0;
    // Если изолированный объект меньше этой доли от локального эталонного размера, он считается скорее шумом.
    // Увеличение: одиночные мелкие объекты отсекаются строже.
    // Уменьшение: больше маленьких одиночных объектов проходит. Разумно: 0.40-0.90.
    private static final double TUNE_ISOLATED_SMALL_AREA_SCALE = 0.62;
    // Минимальная площадь соседней поддержки, используемая при проверке, одинокий это шум или объект внутри реального поля колоний.
    // Увеличение: от маленьких объектов требуется более сильное соседнее окружение.
    // Уменьшение: проще пропускать изолированные слабые объекты. Разумно: 0.50-1.00.
    private static final double TUNE_ISOLATED_SUPPORT_AREA_SCALE = 0.72;

    // КЛЮЧЕВОЙ ПАРАМЕТР: ОСНОВНОЙ ПОРОГ РАЗМЕРА ДЛЯ TIFF-СКАНОВ ОТНОСИТЕЛЬНО ЛОКАЛЬНОГО ТИПИЧНОГО РАЗМЕРА КОЛОНИИ.
    // Увеличение: фильтр размера для сканов становится строже.
    // Уменьшение: проходит больше мелких и средних колоний, но и больше шумовых объектов. Разумно: 0.75-1.20.
    private static final double TUNE_MIN_AREA_SCALE_SCAN = 0.76;


    // Такой же относительный порог размера, но только для сфотографированных планшетов, а не TIFF-сканов.
    // Увеличение: фильтр для фото становится строже.
    // Уменьшение: проходит больше мелких объектов на фото. Разумно: 0.35-0.90.
    private static final double TUNE_MIN_AREA_SCALE_PHOTO = 0.55;
    // Какая доля слипшегося объекта должна быть покрыта валидными частями после разбиения, чтобы мы приняли разбиение.
    // Увеличение: разбиение принимается только если оно очень убедительное.
    // Уменьшение: легче принять разбиение, но можно получить лишнее дробление. Разумно: 0.45-0.85.
    private static final double TUNE_SPLIT_COVERAGE_FRACTION = 0.82;
    // Максимальный размер, до которого уменьшается большой скан перед анализом.
    // Увеличение: анализ идет на более детальном изображении, но дольше и тяжелее по памяти.
    // Уменьшение: быстрее работа, но теряется мелкая детализация. Разумно: 1400-5000.
    private static final double TUNE_SCAN_MAX_DIM = 4000.0;
    // Насколько квадратный crop вокруг лунки должен быть больше, чем найденный круг лунки.
    // Увеличение: вокруг лунки захватывается больше запаса.
    // Уменьшение: crop плотнее, но можно срезать краевые колонии. Разумно: 1.05-1.40.
    private static final double TUNE_WELL_CROP_SCALE = 1.28;
    // Горизонтальное положение ROI скана внутри полного изображения; помогает быстрее выделить область планшета.
    // Увеличение: ROI смещается правее.
    // Уменьшение: ROI смещается левее. Держать в пределах 0.0-1.0, обычно 0.40-0.65.
    private static final double TUNE_SCAN_ROI_LEFT_FRACTION = 0.52;
    // Вертикальное положение ROI скана внутри полного изображения; помогает быстрее выделить область планшета.
    // Увеличение: ROI смещается ниже.
    // Уменьшение: ROI смещается выше. Держать в пределах 0.0-1.0, обычно 0.15-0.40.
    private static final double TUNE_SCAN_ROI_TOP_FRACTION = 0.26;

    // КЛЮЧЕВОЙ ПАРАМЕТР: МНОЖИТЕЛЬ ПОРОГА SCORE ДЛЯ TIFF-СКАНОВ.
    // Увеличение: сегментация сканов строже, слабые объекты исчезают.
    // Уменьшение: сегментация мягче, появляется больше кандидатов и больше риска мусора. Разумно: 0.85-1.10.
    private static final double TUNE_SCAN_STRONG_THRESHOLD_FACTOR = 0.95;
    // ДОПОЛНИТЕЛЬНЫЙ МЯГКИЙ ПОРОГ ДЛЯ ПОИСКА КРУПНЫХ ПРОПУЩЕННЫХ КОЛОНИЙ НА СКАНАХ.
    // Увеличение: rescue-pass строже и почти ничего не добавляет.
    // Уменьшение: rescue-pass подхватывает больше слабых крупных объектов, но может тащить мусор. Разумно: 0.75-0.92.
    private static final double TUNE_SCAN_RESCUE_THRESHOLD_FACTOR = 0.83;
    // Насколько большой должен быть rescue-кандидат относительно локального типичного размера колонии.
    // Увеличение: rescue-pass берет только очень крупные объекты.
    // Уменьшение: будет добавляться больше средних объектов. Разумно: 0.70-1.10.
    private static final double TUNE_SCAN_RESCUE_MIN_AREA_SCALE = 0.78;
    // Минимальная компактность rescue-кандидата; помогает не тащить рваные и странные пятна.
    // Увеличение: rescue-pass строже к форме.
    // Уменьшение: проходит больше неровных объектов. Разумно: 0.60-0.90.
    private static final double TUNE_SCAN_RESCUE_MIN_SOLIDITY = 0.72;
    // Минимальная округлость rescue-кандидата.
    // Увеличение: проходят в основном более круглые кандидаты.
    // Уменьшение: можно подхватывать более неровные, но отдельные колонии. Разумно: 0.30-0.75.
    private static final double TUNE_SCAN_RESCUE_MIN_CIRCULARITY = 0.38;
    // Максимальная доля перекрытия rescue-кандидата с уже найденными колониями.
    // Увеличение: rescue-pass легче добавляет кандидаты рядом с уже найденными.
    // Уменьшение: rescue-pass берет почти только полностью новые объекты. Разумно: 0.05-0.35.
    private static final double TUNE_SCAN_RESCUE_MAX_OVERLAP = 0.18;
    // Минимальная доля новых пикселей у rescue-кандидата относительно уже найденных объектов.
    // Увеличение: берутся только действительно новые объекты.
    // Уменьшение: можно добавлять кандидаты, частично сливающиеся с уже найденными. Разумно: 0.35-0.85.
    private static final double TUNE_SCAN_RESCUE_MIN_NEW_FRACTION = 0.55;
    // Минимальный средний score rescue-кандидата, если форма у него неидеальная.
    // Увеличение: rescue-pass строже к слабоконтрастным объектам.
    // Уменьшение: легче подхватываются бледные, но большие колонии. Разумно: 0.04-0.10.
    private static final double TUNE_SCAN_RESCUE_MIN_SCORE = 0.058;
    // Максимальная средняя яркость rescue-кандидата; слишком светлые объекты считаются менее надежными.
    // Увеличение: rescue-pass допускает более светлые кандидаты.
    // Уменьшение: спасаются только более темные объекты. Разумно: 0.14-0.24.
    private static final double TUNE_SCAN_RESCUE_MAX_GRAY = 0.20;
    // Минимальная фиолетовая окраска rescue-кандидата.
    // Увеличение: требуется более выраженная окраска.
    // Уменьшение: rescue-pass допускает более бледные объекты. Разумно: 0.015-0.07.
    private static final double TUNE_SCAN_RESCUE_MIN_PURPLE = 0.02;
    // Минимальная заполненность слабой маски, чтобы несколько уже найденных сегментов можно было склеить обратно в один blob.
    // Увеличение: merge-back будет срабатывать только на очень цельных пятнах.
    // Уменьшение: merge-back начнет чаще склеивать соседние сегменты. Разумно: 0.75-0.95.
    private static final double TUNE_SCAN_MERGE_MIN_SOLIDITY = 0.84;
    // Минимальная округлость слабой маски для merge-back.
    // Увеличение: склеиваются в основном круглые пятна.
    // Уменьшение: можно склеивать и более неровные пятна. Разумно: 0.20-0.70.
    private static final double TUNE_SCAN_MERGE_MIN_CIRCULARITY = 0.30;
    // Максимальная вытянутость слабой маски для merge-back.
    // Увеличение: вытянутые пятна тоже могут быть склеены.
    // Уменьшение: merge-back работает только на более компактных blob. Разумно: 1.2-2.0.
    private static final double TUNE_SCAN_MERGE_MAX_ASPECT = 1.65;
    // Сколько уже найденных сегментов минимум должна накрывать слабая маска, чтобы их имело смысл схлопывать в один объект.
    // Увеличение: merge-back будет реже.
    // Уменьшение: можно склеивать даже пары сегментов. Разумно: 2-5.
    private static final int TUNE_SCAN_MERGE_MIN_SEGMENTS = 3;
    // Какая доля площади слабой маски должна быть уже покрыта найденными сегментами, чтобы считать это одним blob, а не случайным объединением.
    // Увеличение: merge-back строже.
    // Уменьшение: merge-back будет срабатывать чаще. Разумно: 0.25-0.75.
    private static final double TUNE_SCAN_MERGE_MIN_COVERAGE = 0.40;

    // Нижний квантиль, который используется при оценке "типичного" размера колонии.
    // Увеличение: типичный размер сдвигается вверх, и фильтр размера становится строже.
    // Уменьшение: типичный размер сдвигается вниз, и больше средних колоний проходит. Разумно: 0.10-0.35.
    private static final double TUNE_SINGLE_AREA_LOW_QUANTILE = 0.18;
    // Верхний квантиль, который используется при оценке "типичного" размера колонии.
    // Увеличение: в расчет сильнее попадают крупные и слипшиеся объекты.
    // Уменьшение: влияние крупных объектов уменьшается. Разумно: 0.55-0.85.
    private static final double TUNE_SINGLE_AREA_HIGH_QUANTILE = 0.68;
    // Какой квантиль внутри очищенного диапазона считать итоговым "типичным" размером колонии.
    // Увеличение: типичный размер больше, фильтр размера жестче.
    // Уменьшение: типичный размер меньше, проходит больше колоний. Разумно: 0.20-0.55.
    private static final double TUNE_SINGLE_AREA_TARGET_QUANTILE = 0.35;

    // Насколько большим должен быть объект относительно типичного размера, чтобы мы попытались его разделить.
    // Увеличение: крупные слипшиеся пятна разделяются реже.
    // Уменьшение: разделение пробуется чаще. Разумно: 1.30-2.20.
    private static final double TUNE_SPLIT_TRIGGER_AREA_SCALE = 1.35;
    // Минимальный размер дочерней части после разделения относительно типичного размера колонии.
    // Увеличение: мелкие части после разделения чаще отбрасываются.
    // Уменьшение: разделение принимает более мелкие части. Разумно: 0.45-0.90.
    private static final double TUNE_SPLIT_MIN_CHILD_AREA_SCALE = 0.45;
    // Минимальное покрытие родительского объекта валидными дочерними частями.
    // Увеличение: разделение принимается только если оно очень убедительное.
    // Уменьшение: разделение проходит легче. Разумно: 0.35-0.75.
    private static final double TUNE_SPLIT_VALID_COVERAGE = 0.35;
    // Если родительское пятно достаточно цельное и компактное, то split включается только для очень больших объектов.
    // Увеличение: цельные большие пятна реже дробятся.
    // Уменьшение: даже компактные большие пятна будут дробиться чаще. Разумно: 10-30.
    private static final double TUNE_SPLIT_COMPACT_TRIGGER_AREA_SCALE = 18.0;
    // Минимальная заполненность (solidity), чтобы считать пятно компактным и не спешить его дробить.
    // Увеличение: защита от дробления действует только на очень цельные пятна.
    // Уменьшение: защита будет срабатывать чаще. Разумно: 0.80-0.98.
    private static final double TUNE_SPLIT_COMPACT_MIN_SOLIDITY = 0.90;
    // Минимальная округлость, чтобы считать пятно компактным и, скорее всего, одиночной колонией.
    // Увеличение: защита действует только на более круглые пятна.
    // Уменьшение: защита от дробления применяется шире. Разумно: 0.35-0.80.
    private static final double TUNE_SPLIT_COMPACT_MIN_CIRCULARITY = 0.48;
    // Максимальная вытянутость компактного пятна по bounding box.
    // Увеличение: вытянутые пятна тоже реже будут дробиться.
    // Уменьшение: защита останется только для более круглых пятен. Разумно: 1.2-2.2.
    private static final double TUNE_SPLIT_COMPACT_MAX_ASPECT = 1.55;
    // Насколько далеко должны разойтись центры дочерних частей, чтобы split считался правдоподобным.
    // Увеличение: требуется более явное разделение частей.
    // Уменьшение: split принимается даже для очень близких частей. Разумно: 0.7-1.8.
    private static final double TUNE_SPLIT_MIN_CENTROID_DISTANCE_SCALE = 1.05;

    // Маленькие слабые объекты меньше этого множителя от абсолютного минимума считаются вероятным мусором.
    // Увеличение: правило отсекает больше маленьких слабых объектов.
    // Уменьшение: больше мелочи проходит дальше. Разумно: 2.0-4.5.
    private static final double TUNE_SCAN_SMALL_WEAK_AREA_MULTIPLIER = 3.8;
    // Маленькие слабые объекты меньше этой доли от локального размера колонии считаются вероятным мусором.
    // Увеличение: фильтр мелких слабых объектов строже.
    // Уменьшение: мелкие кандидаты чаще выживают. Разумно: 0.70-1.20.
    private static final double TUNE_SCAN_SMALL_WEAK_LOCAL_SCALE = 0.90;
    // Маленькие слабые объекты допускаются только если у них есть хотя бы такая поддержка соседями.
    // Увеличение: одиночным слабым точкам сложнее пройти.
    // Уменьшение: слабые объекты легче проходят даже без соседей. Разумно: 0-3.
    private static final int TUNE_SCAN_SMALL_WEAK_MAX_SUPPORT = 1;
    // Если маленький объект светлее этого порога, он выглядит подозрительно и чаще отбрасывается.
    // Увеличение: фильтр по яркости слабее, светлые объекты чаще проходят.
    // Уменьшение: фильтр по яркости строже. Разумно: 0.12-0.28.
    private static final double TUNE_SCAN_SMALL_WEAK_MAX_GRAY = 0.2;
    // Если у маленького объекта меньше фиолетового сигнала, чем здесь, он меньше похож на окрашенную колонию.
    // Увеличение: требуется более сильная окраска, фильтр строже.
    // Уменьшение: допускаются более бледные объекты. Разумно: 0.02-0.08.
    private static final double TUNE_SCAN_SMALL_WEAK_MIN_PURPLE = 0.02;
    // Минимальный внутренний score, при котором маленький слабый объект все еще может быть оставлен.
    // Увеличение: спасается меньше мелких спорных объектов.
    // Уменьшение: спасается больше мелких спорных объектов. Разумно: 0.04-0.12.
    private static final double TUNE_SCAN_SMALL_WEAK_MIN_SCORE = 0.075;
    // Множитель локального порога в правиле отсева маленьких слабых объектов.
    // Увеличение: маленьким слабым объектам сложнее пройти.
    // Уменьшение: правило мягче. Разумно: 1.5-3.0.
    private static final double TUNE_SCAN_SMALL_WEAK_THRESHOLD_SCALE = 1.35;
    // Пограничный объект все еще может быть принят, если он не меньше этой доли от локально ожидаемого размера колонии.
    // Увеличение: rescue для пограничных объектов строже по размеру.
    // Уменьшение: rescue срабатывает чаще. Разумно: 0.60-0.95.
    private static final double TUNE_BORDERLINE_AREA_SCALE = 0.72;
    // Пограничные объекты должны быть достаточно компактными и заполненными, а не рваными.
    // Увеличение: rescue строже к форме объекта.
    // Уменьшение: можно пропускать более рыхлые объекты. Разумно: 0.75-0.98.
    private static final double TUNE_BORDERLINE_MIN_SOLIDITY = 0.8;
    // Пограничные объекты не должны быть слишком вытянутыми; уменьшить — чтобы принимать более неровные формы.
    // Увеличение: принимаются только более круглые объекты.
    // Уменьшение: принимаются более вытянутые и неровные формы. Разумно: 0.45-0.90.
    private static final double TUNE_BORDERLINE_MIN_CIRCULARITY = 0.52;
    // Пограничным объектам проще доверять, если они достаточно темные.
    // Увеличение: rescue допускает более светлые объекты.
    // Уменьшение: rescue строже, требуется более темный объект. Разумно: 0.10-0.25.
    private static final double TUNE_BORDERLINE_MAX_GRAY = 0.20;
    // Пограничным объектам проще доверять, если у них достаточно фиолетовой окраски.
    // Увеличение: требуется более сильная окраска.
    // Уменьшение: rescue допускает более бледные объекты. Разумно: 0.02-0.08.
    private static final double TUNE_BORDERLINE_MIN_PURPLE = 0.05;

    // КЛЮЧЕВОЙ ПАРАМЕТР: МИНИМАЛЬНЫЙ ВНУТРЕННИЙ SCORE ДЛЯ "СПАСЕНИЯ" ПОГРАНИЧНОГО ОБЪЕКТА.
    // Увеличение: спасается меньше сомнительных колоний.
    // Уменьшение: rescue срабатывает чаще и может вернуть спорные объекты. Разумно: 0.06-0.14.
    private static final double TUNE_BORDERLINE_MIN_SCORE = 0.082;
    
    // Множитель локального порога в правиле "спасения" пограничных объектов.
    // Увеличение: rescue становится строже.
    // Уменьшение: rescue становится мягче. Разумно: 1.8-3.5.
    private static final double TUNE_BORDERLINE_THRESHOLD_SCALE = 2.7;

    private static Config CONFIG = new Config();
    private static List<LayoutSeed> CUSTOM_LAYOUT = List.of();
    private static ManifestIndex MANIFEST = ManifestIndex.empty();
    private static List<String> SUMMARY_GROUP_FIELDS = List.of();
    private static DatasetMode EFFECTIVE_DATASET_MODE = DatasetMode.AUTO;

    private record Config(
        double backgroundSigmaScale,
        double scoreQuantile,
        int minColonyArea,
        double thresholdRelaxation,
        double clumpAreaFactor,
        double neighborRadiusScale,
        int countMaskErosionPx,
        double outerRejectBandPx,
        double isolatedSmallAreaScale,
        double isolatedSupportAreaScale,
        double minAreaScaleScan,
        double minAreaScalePhoto,
        double splitCoverageFraction,
        double scanMaxDim,
        double wellCropScale,
        double scanRoiLeftFraction,
        double scanRoiTopFraction,
        double scanThresholdFactor
    ) {
        Config() {
            this(
                TUNE_BACKGROUND_SIGMA_SCALE,
                TUNE_SCORE_QUANTILE,
                TUNE_MIN_COLONY_AREA,
                TUNE_THRESHOLD_RELAXATION,
                TUNE_CLUMP_AREA_FACTOR,
                TUNE_NEIGHBOR_RADIUS_SCALE,
                TUNE_COUNT_MASK_EROSION_PX,
                TUNE_OUTER_REJECT_BAND_PX,
                TUNE_ISOLATED_SMALL_AREA_SCALE,
                TUNE_ISOLATED_SUPPORT_AREA_SCALE,
                TUNE_MIN_AREA_SCALE_SCAN,
                TUNE_MIN_AREA_SCALE_PHOTO,
                TUNE_SPLIT_COVERAGE_FRACTION,
                TUNE_SCAN_MAX_DIM,
                TUNE_WELL_CROP_SCALE,
                TUNE_SCAN_ROI_LEFT_FRACTION,
                TUNE_SCAN_ROI_TOP_FRACTION,
                TUNE_SCAN_STRONG_THRESHOLD_FACTOR
            );
        }

        static Config fromProperties(Properties properties) {
            Config defaults = new Config();
            return new Config(
                getDouble(properties, "background_sigma_scale", defaults.backgroundSigmaScale()),
                getDouble(properties, "score_quantile", defaults.scoreQuantile()),
                getInt(properties, "min_colony_area", defaults.minColonyArea()),
                getDouble(properties, "threshold_relaxation", defaults.thresholdRelaxation()),
                getDouble(properties, "clump_area_factor", defaults.clumpAreaFactor()),
                getDouble(properties, "neighbor_radius_scale", defaults.neighborRadiusScale()),
                getInt(properties, "count_mask_erosion_px", defaults.countMaskErosionPx()),
                getDouble(properties, "outer_reject_band_px", defaults.outerRejectBandPx()),
                getDouble(properties, "isolated_small_area_scale", defaults.isolatedSmallAreaScale()),
                getDouble(properties, "isolated_support_area_scale", defaults.isolatedSupportAreaScale()),
                getDouble(properties, "min_area_scale_scan", defaults.minAreaScaleScan()),
                getDouble(properties, "min_area_scale_photo", defaults.minAreaScalePhoto()),
                getDouble(properties, "split_coverage_fraction", defaults.splitCoverageFraction()),
                getDouble(properties, "scan_max_dim", defaults.scanMaxDim()),
                getDouble(properties, "well_crop_scale", defaults.wellCropScale()),
                getDouble(properties, "scan_roi_left_fraction", defaults.scanRoiLeftFraction()),
                getDouble(properties, "scan_roi_top_fraction", defaults.scanRoiTopFraction()),
                getDouble(properties, "scan_threshold_factor", defaults.scanThresholdFactor())
            );
        }
    }

    private enum DatasetMode {
        AUTO,
        HF,
        ZR;

        static DatasetMode parse(String raw) {
            String normalized = raw.trim().toLowerCase(Locale.ROOT);
            return switch (normalized) {
                case "auto" -> AUTO;
                case "hf" -> HF;
                case "zr" -> ZR;
                default -> throw new IllegalArgumentException("Unsupported --dataset value: " + raw + ". Expected one of: auto, hf, zr");
            };
        }
    }

    private record Args(
        Path outputDir,
        DatasetMode datasetMode,
        Path settingsPath,
        Path layoutPath,
        Path manifestPath,
        List<Path> inputs
    ) {}

    private record Circle(double x, double y, double radius, double score) {}

    private record LayoutSeed(int wellIndex, double xFraction, double yFraction, double radiusFraction) {}

    private record Box(int left, int top, int right, int bottom) {
        int width() {
            return right - left;
        }

        int height() {
            return bottom - top;
        }
    }

    private record Metadata(
        String imageName,
        String plateKey,
        int captureIndex,
        Double doseGy,
        String doseUnit,
        Double topConcentration,
        Double bottomConcentration,
        boolean reverseColumnAssignments,
        Map<String, Object> imageFields,
        Map<Integer, Map<String, Object>> wellFields
    ) {}

    private record WellAssignment(
        String rowName,
        Double concentration,
        int replicate,
        Map<String, Object> extraFields
    ) {}

    private record ManifestIndex(
        Map<String, Map<Integer, Map<String, Object>>> byImageName,
        Map<String, Map<Integer, Map<String, Object>>> byImageStem
    ) {
        static ManifestIndex empty() {
            return new ManifestIndex(Map.of(), Map.of());
        }

        Map<Integer, Map<String, Object>> rowsForImage(String imageName) {
            Map<Integer, Map<String, Object>> direct = byImageName.get(imageName);
            if (direct != null) {
                return direct;
            }
            return byImageStem.getOrDefault(stripExtension(imageName), Map.of());
        }
    }

    private record AnalysisFrame(
        BufferedImage image,
        Path imagePath,
        String analysisSource,
        Box analysisRoiBox,
        Box analysisCanvas,
        double analysisScale
    ) {}

    private record WellRegion(
        int wellIndex,
        Circle roughCircle,
        Circle refinedCircle,
        Box cropBox,
        boolean[][] maskLocal,
        boolean[][] coreMaskLocal
    ) {}

    private record RawParticle(
        int id,
        Roi roiLocal,
        int area,
        double perimeter,
        double circularity,
        double solidity,
        double eccentricity,
        double centroidX,
        double centroidY,
        double meanScore,
        double meanGray,
        double meanPurple,
        double edgeDistance,
        double coreOverlap,
        Rectangle bounds
    ) {
        double bboxAspectRatio() {
            double width = Math.max(1.0, bounds.width);
            double height = Math.max(1.0, bounds.height);
            return Math.max(width, height) / Math.min(width, height);
        }

        double bboxExtent() {
            double bboxArea = Math.max(1.0, bounds.width * bounds.height);
            return area / bboxArea;
        }
    }

    private record SegmentRecord(
        int wellIndex,
        String rowName,
        Double concentration,
        int replicate,
        Map<String, Object> extraFields,
        int areaPx,
        double eccentricity,
        double circularity,
        double solidity,
        double edgeDistancePx,
        double coreOverlap,
        double meanGray,
        double meanPurple,
        double centroidX,
        double centroidY,
        Rectangle bbox,
        Roi roiGlobal,
        int labelIndex
    ) {}

    private record WellSummary(
        int wellIndex,
        String rowName,
        int replicate,
        Double concentration,
        Map<String, Object> extraFields,
        int rawComponentCount,
        int countedColonies,
        double singleAreaPx,
        double threshold,
        int wellAreaPx,
        int countAreaPx
    ) {}

    private record WellOutput(
        List<SegmentRecord> segments,
        WellSummary summary,
        boolean[][] countedMaskLocal
    ) {}

    private record ConnectedComponent(List<int[]> pixels, Rectangle bounds) {}

    private record ComponentShape(
        double circularity,
        double solidity,
        double aspectRatio,
        double centroidX,
        double centroidY
    ) {}

    public static void main(String[] args) throws Exception {
        Args parsed = parseArgs(args);
        applyExternalConfiguration(parsed);
        Files.createDirectories(parsed.outputDir);

        List<Map<String, Object>> allWellRows = new ArrayList<>();
        List<Map<String, Object>> runSummaries = new ArrayList<>();

        for (Path input : parsed.inputs) {
            Map<String, Object> summary = analyzeImage(input, parsed.outputDir, EFFECTIVE_DATASET_MODE, allWellRows);
            runSummaries.add(summary);
        }

        writeDuplicateWellSummary(allWellRows, parsed.outputDir.resolve("duplicate_well_means.csv"));
        writeDuplicateGroupSummary(allWellRows, parsed.outputDir.resolve("duplicate_group_means.csv"));

        int totalCount = 0;
        for (Map<String, Object> run : runSummaries) {
            totalCount += ((Number) run.get("total_counted_colonies")).intValue();
        }

        Path runSummaryPath = parsed.outputDir.resolve("run_summary.json");
        try (BufferedWriter writer = Files.newBufferedWriter(runSummaryPath, StandardCharsets.UTF_8)) {
            writer.write("{\n");
            writer.write("  \"runs\": [\n");
            for (int i = 0; i < runSummaries.size(); i++) {
                writer.write(toJson(runSummaries.get(i), 4));
                if (i + 1 < runSummaries.size()) {
                    writer.write(",");
                }
                writer.write("\n");
            }
            writer.write("  ],\n");
            writer.write("  \"totals\": {\n");
            writer.write("    \"images_processed\": " + runSummaries.size() + ",\n");
            writer.write("    \"total_counted_colonies\": " + totalCount + "\n");
            writer.write("  }\n");
            writer.write("}\n");
        }
    }

    private static Args parseArgs(String[] args) throws IOException {
        Path outputDir = Path.of("analysis_output_java");
        DatasetMode datasetMode = DatasetMode.AUTO;
        Path settingsPath = null;
        Path layoutPath = null;
        Path manifestPath = null;
        List<String> rawInputs = new ArrayList<>();
        for (int i = 0; i < args.length; i++) {
            if ("--output-dir".equals(args[i]) && i + 1 < args.length) {
                outputDir = Path.of(args[++i]).toAbsolutePath().normalize();
            } else if ("--dataset".equals(args[i]) && i + 1 < args.length) {
                datasetMode = DatasetMode.parse(args[++i]);
            } else if ("--config".equals(args[i]) && i + 1 < args.length) {
                settingsPath = Path.of(args[++i]).toAbsolutePath().normalize();
            } else if ("--layout".equals(args[i]) && i + 1 < args.length) {
                layoutPath = Path.of(args[++i]).toAbsolutePath().normalize();
            } else if ("--manifest".equals(args[i]) && i + 1 < args.length) {
                manifestPath = Path.of(args[++i]).toAbsolutePath().normalize();
            } else {
                rawInputs.add(args[i]);
            }
        }
        if (rawInputs.isEmpty()) {
            throw new IllegalArgumentException("Usage: java ClonogenicAnalyzer [--output-dir DIR] [--dataset auto|hf|zr] [--config settings.properties] [--layout well_layout.csv] [--manifest plate_manifest.csv] <images...>");
        }

        List<Path> inputs = new ArrayList<>();
        for (String raw : rawInputs) {
            if (raw.contains("*") || raw.contains("?")) {
                Path parent = Path.of(raw).toAbsolutePath().normalize().getParent();
                if (parent == null) {
                    parent = Paths.get("").toAbsolutePath().normalize();
                }
                String pattern = Path.of(raw).getFileName().toString();
                try (DirectoryStream<Path> stream = Files.newDirectoryStream(parent, pattern)) {
                    for (Path path : stream) {
                        inputs.add(path.toAbsolutePath().normalize());
                    }
                }
            } else {
                inputs.add(Path.of(raw).toAbsolutePath().normalize());
            }
        }
        inputs.sort(Comparator.naturalOrder());
        if (inputs.isEmpty()) {
            throw new IllegalArgumentException("No input images matched the provided paths.");
        }
        return new Args(outputDir, datasetMode, settingsPath, layoutPath, manifestPath, inputs);
    }

    private static void applyExternalConfiguration(Args parsed) throws IOException {
        CONFIG = new Config();
        CUSTOM_LAYOUT = List.of();
        MANIFEST = ManifestIndex.empty();
        SUMMARY_GROUP_FIELDS = List.of();
        EFFECTIVE_DATASET_MODE = parsed.datasetMode();

        if (parsed.settingsPath() != null) {
            Properties properties = new Properties();
            try (BufferedReader reader = Files.newBufferedReader(parsed.settingsPath(), StandardCharsets.UTF_8)) {
                properties.load(reader);
            }
            CONFIG = Config.fromProperties(properties);
            SUMMARY_GROUP_FIELDS = parseListProperty(properties.getProperty("summary_group_fields"));
            if (EFFECTIVE_DATASET_MODE == DatasetMode.AUTO) {
                String configuredMode = properties.getProperty("dataset_mode");
                if (configuredMode != null && !configuredMode.isBlank()) {
                    EFFECTIVE_DATASET_MODE = DatasetMode.parse(configuredMode);
                }
            }
        }

        if (parsed.layoutPath() != null) {
            CUSTOM_LAYOUT = loadLayoutSeeds(parsed.layoutPath());
        }
        if (parsed.manifestPath() != null) {
            MANIFEST = loadManifest(parsed.manifestPath());
        }
    }

    private static List<String> parseListProperty(String raw) {
        if (raw == null || raw.isBlank()) {
            return List.of();
        }
        String[] parts = raw.split(",");
        List<String> values = new ArrayList<>();
        for (String part : parts) {
            String trimmed = part.trim();
            if (!trimmed.isEmpty()) {
                values.add(trimmed);
            }
        }
        return values;
    }

    private static double getDouble(Properties properties, String key, double defaultValue) {
        String raw = properties.getProperty(key);
        if (raw == null || raw.isBlank()) {
            return defaultValue;
        }
        return Double.parseDouble(raw.trim().replace(',', '.'));
    }

    private static int getInt(Properties properties, String key, int defaultValue) {
        String raw = properties.getProperty(key);
        if (raw == null || raw.isBlank()) {
            return defaultValue;
        }
        return Integer.parseInt(raw.trim());
    }

    private static List<LayoutSeed> loadLayoutSeeds(Path path) throws IOException {
        CsvTable table = readCsvTable(path);
        List<LayoutSeed> seeds = new ArrayList<>();
        for (Map<String, String> row : table.rows()) {
            int wellIndex = parseRequiredInt(row, "well_index", path);
            double xFraction = parseRequiredDouble(row, "x_fraction", path);
            double yFraction = parseRequiredDouble(row, "y_fraction", path);
            double radiusFraction = parseRequiredDouble(row, "radius_fraction", path);
            seeds.add(new LayoutSeed(wellIndex, xFraction, yFraction, radiusFraction));
        }
        seeds.sort(Comparator.comparingInt(LayoutSeed::wellIndex));
        return seeds;
    }

    private static ManifestIndex loadManifest(Path path) throws IOException {
        CsvTable table = readCsvTable(path);
        Map<String, Map<Integer, Map<String, Object>>> byImageName = new LinkedHashMap<>();
        Map<String, Map<Integer, Map<String, Object>>> byImageStem = new LinkedHashMap<>();
        for (Map<String, String> row : table.rows()) {
            String imageName = row.getOrDefault("image_name", "").trim();
            if (imageName.isEmpty()) {
                throw new IllegalArgumentException("Manifest row without image_name in " + path);
            }
            int wellIndex = parseRequiredInt(row, "well_index", path);
            Map<String, Object> values = new LinkedHashMap<>();
            for (Map.Entry<String, String> entry : row.entrySet()) {
                String key = entry.getKey();
                String value = entry.getValue();
                if (value == null || value.isBlank()) {
                    continue;
                }
                values.put(key, value.trim());
            }
            byImageName.computeIfAbsent(imageName, unused -> new LinkedHashMap<>()).put(wellIndex, values);
            byImageStem.computeIfAbsent(stripExtension(imageName), unused -> new LinkedHashMap<>()).put(wellIndex, values);
        }
        return new ManifestIndex(byImageName, byImageStem);
    }

    private record CsvTable(List<String> headers, List<Map<String, String>> rows) {}

    private static CsvTable readCsvTable(Path path) throws IOException {
        List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
        if (lines.isEmpty()) {
            throw new IllegalArgumentException("Empty CSV file: " + path);
        }
        char delimiter = detectDelimiter(lines.getFirst());
        List<String> headers = parseCsvLine(lines.getFirst(), delimiter);
        List<Map<String, String>> rows = new ArrayList<>();
        for (int i = 1; i < lines.size(); i++) {
            String line = lines.get(i);
            if (line == null || line.isBlank()) {
                continue;
            }
            List<String> values = parseCsvLine(line, delimiter);
            Map<String, String> row = new LinkedHashMap<>();
            for (int col = 0; col < headers.size(); col++) {
                String value = col < values.size() ? values.get(col) : "";
                row.put(headers.get(col).trim(), value.trim());
            }
            rows.add(row);
        }
        return new CsvTable(headers, rows);
    }

    private static char detectDelimiter(String header) {
        long semicolons = header.chars().filter(ch -> ch == ';').count();
        long commas = header.chars().filter(ch -> ch == ',').count();
        return semicolons > commas ? ';' : ',';
    }

    private static List<String> parseCsvLine(String line, char delimiter) {
        List<String> out = new ArrayList<>();
        StringBuilder current = new StringBuilder();
        boolean inQuotes = false;
        for (int i = 0; i < line.length(); i++) {
            char ch = line.charAt(i);
            if (ch == '"') {
                if (inQuotes && i + 1 < line.length() && line.charAt(i + 1) == '"') {
                    current.append('"');
                    i++;
                } else {
                    inQuotes = !inQuotes;
                }
            } else if (ch == delimiter && !inQuotes) {
                out.add(current.toString());
                current.setLength(0);
            } else {
                current.append(ch);
            }
        }
        out.add(current.toString());
        return out;
    }

    private static int parseRequiredInt(Map<String, String> row, String key, Path path) {
        String value = row.getOrDefault(key, "").trim();
        if (value.isEmpty()) {
            throw new IllegalArgumentException("Missing required integer field '" + key + "' in " + path);
        }
        return Integer.parseInt(value);
    }

    private static double parseRequiredDouble(Map<String, String> row, String key, Path path) {
        String value = row.getOrDefault(key, "").trim();
        if (value.isEmpty()) {
            throw new IllegalArgumentException("Missing required numeric field '" + key + "' in " + path);
        }
        return Double.parseDouble(value.replace(',', '.'));
    }

    private static Map<String, Object> analyzeImage(
        Path imagePath,
        Path outputRoot,
        DatasetMode datasetMode,
        List<Map<String, Object>> allWellRows
    ) throws Exception {
        Path imageOutputDir = outputRoot.resolve(stripExtension(imagePath.getFileName().toString()));
        Files.createDirectories(imageOutputDir);

        Metadata metadata = parseMetadata(imagePath, datasetMode);
        BufferedImage source = openBuffered(imagePath);
        AnalysisFrame frame = prepareAnalysisFrame(imagePath, source, imageOutputDir);

        double[][] gray = toGray(frame.image());
        double[][] purple = toPurple(frame.image());
        List<Circle> roughSeeds;
        if (frame.analysisSource().startsWith("scan_roi")) {
            roughSeeds = createScanSeedCircles(frame.image().getWidth(), frame.image().getHeight());
        } else {
            Box plateBox = detectPlateBox(gray);
            roughSeeds = createSeedCircles(plateBox, frame.image().getWidth(), frame.image().getHeight());
        }
        writeSeedCsv(imageOutputDir.resolve("java_seed_circles.csv"), roughSeeds);
        List<Circle> refined = refineWells(frame.imagePath(), imageOutputDir.resolve("java_seed_circles.csv"), imageOutputDir.resolve("java_refined_wells.csv"));

        List<WellRegion> wellRegions = buildWellRegions(refined, frame.image().getWidth(), frame.image().getHeight());
        double meanRadius = refined.stream().mapToDouble(Circle::radius).average().orElse(80.0);
        double[][] scoreMap = computeScoreMap(gray, purple, meanRadius);

        List<SegmentRecord> allSegments = new ArrayList<>();
        List<WellSummary> wellSummaries = new ArrayList<>();
        boolean[][] globalMask = new boolean[frame.image().getHeight()][frame.image().getWidth()];

        for (WellRegion region : wellRegions) {
            WellOutput output = analyzeWell(scoreMap, gray, purple, region, metadata, frame.analysisSource().startsWith("scan_roi"));
            allSegments.addAll(output.segments());
            wellSummaries.add(output.summary());
            mergeLocalMask(globalMask, output.countedMaskLocal(), region.cropBox());
            saveWellDebug(frame.image(), scoreMap, region, output, imageOutputDir);
        }

        savePlateOverlay(frame.image(), wellRegions, allSegments, wellSummaries, imageOutputDir.resolve("plate_overlay.png"));
        saveColoniesOnly(frame.image(), globalMask, imageOutputDir.resolve("colonies_only.png"));

        List<Map<String, Object>> wellRows = new ArrayList<>();
        for (int i = 0; i < wellRegions.size(); i++) {
            WellRegion region = wellRegions.get(i);
            WellSummary summary = wellSummaries.get(i);
            Map<String, Object> row = baseMetadataMap(metadata);
            row.putAll(summary.extraFields());
            row.put("well_index", summary.wellIndex());
            row.put("row_name", summary.rowName());
            row.put("replicate", summary.replicate());
            row.put("concentration", summary.concentration());
            row.put("raw_component_count", summary.rawComponentCount());
            row.put("counted_colonies", summary.countedColonies());
            row.put("single_area_px", round(summary.singleAreaPx(), 3));
            row.put("threshold", round(summary.threshold(), 6));
            row.put("well_area_px", summary.wellAreaPx());
            row.put("count_area_px", summary.countAreaPx());
            row.put("analysis_source", frame.analysisSource());
            row.put("analysis_offset_x_px", frame.analysisCanvas().left());
            row.put("analysis_offset_y_px", frame.analysisCanvas().top());
            row.put("analysis_scale", round(frame.analysisScale(), 6));
            row.put("well_center_x_px", round(region.refinedCircle().x(), 3));
            row.put("well_center_y_px", round(region.refinedCircle().y(), 3));
            row.put("global_well_center_x_px", round(region.refinedCircle().x() / frame.analysisScale() + frame.analysisCanvas().left(), 3));
            row.put("global_well_center_y_px", round(region.refinedCircle().y() / frame.analysisScale() + frame.analysisCanvas().top(), 3));
            row.put("well_radius_px", round(region.refinedCircle().radius(), 3));
            row.put("rough_center_x_px", round(region.roughCircle().x(), 3));
            row.put("rough_center_y_px", round(region.roughCircle().y(), 3));
            row.put("global_rough_center_x_px", round(region.roughCircle().x() / frame.analysisScale() + frame.analysisCanvas().left(), 3));
            row.put("global_rough_center_y_px", round(region.roughCircle().y() / frame.analysisScale() + frame.analysisCanvas().top(), 3));
            row.put("rough_radius_px", round(region.roughCircle().radius(), 3));
            wellRows.add(row);
            allWellRows.add(row);
        }

        List<Map<String, Object>> segmentRows = new ArrayList<>();
        for (int i = 0; i < allSegments.size(); i++) {
            SegmentRecord segment = allSegments.get(i);
            Map<String, Object> row = baseMetadataMap(metadata);
            row.putAll(segment.extraFields());
            row.put("segment_index", i + 1);
            row.put("well_index", segment.wellIndex());
            row.put("row_name", segment.rowName());
            row.put("concentration", segment.concentration());
            row.put("replicate", segment.replicate());
            row.put("label_index", segment.labelIndex());
            row.put("area_px", segment.areaPx());
            row.put("eccentricity", round(segment.eccentricity(), 4));
            row.put("circularity", round(segment.circularity(), 4));
            row.put("solidity", round(segment.solidity(), 4));
            row.put("edge_distance_px", round(segment.edgeDistancePx(), 3));
            row.put("core_overlap", round(segment.coreOverlap(), 4));
            row.put("mean_gray", round(segment.meanGray(), 5));
            row.put("mean_purple", round(segment.meanPurple(), 5));
            row.put("centroid_x_px", round(segment.centroidX(), 3));
            row.put("centroid_y_px", round(segment.centroidY(), 3));
            row.put("global_centroid_x_px", round(segment.centroidX() / frame.analysisScale() + frame.analysisCanvas().left(), 3));
            row.put("global_centroid_y_px", round(segment.centroidY() / frame.analysisScale() + frame.analysisCanvas().top(), 3));
            row.put("bbox_left_px", segment.bbox().x);
            row.put("bbox_top_px", segment.bbox().y);
            row.put("bbox_right_px", segment.bbox().x + segment.bbox().width);
            row.put("bbox_bottom_px", segment.bbox().y + segment.bbox().height);
            segmentRows.add(row);
        }

        writeCsv(imageOutputDir.resolve("well_counts.csv"), wellRows);
        writeCsv(imageOutputDir.resolve("segments.csv"), segmentRows);
        writeGroupSummary(wellRows, imageOutputDir.resolve("group_summary.csv"));

        Map<String, Object> summary = new HashMap<>();
        summary.put("image", metadata.imageName());
        summary.put("output_dir", imageOutputDir.toString());
        summary.put("geometry_source", "java");
        summary.put("analysis_source", frame.analysisSource());
        summary.put("analysis_scale", round(frame.analysisScale(), 6));
        if (frame.analysisRoiBox() != null) {
            Map<String, Object> roi = new HashMap<>();
            roi.put("left", frame.analysisRoiBox().left());
            roi.put("top", frame.analysisRoiBox().top());
            roi.put("right", frame.analysisRoiBox().right());
            roi.put("bottom", frame.analysisRoiBox().bottom());
            summary.put("analysis_roi_box_px", roi);
        } else {
            summary.put("analysis_roi_box_px", null);
        }
        Map<String, Object> canvas = new HashMap<>();
        canvas.put("left", frame.analysisCanvas().left());
        canvas.put("top", frame.analysisCanvas().top());
        summary.put("analysis_canvas_origin_px", canvas);
        summary.put("metadata", metadataToMap(metadata));
        summary.put("well_summaries", wellSummariesToMaps(wellSummaries));
        int totalCount = wellSummaries.stream().mapToInt(WellSummary::countedColonies).sum();
        summary.put("total_counted_colonies", totalCount);
        writeJson(imageOutputDir.resolve("summary.json"), summary);
        return summary;
    }

    private static BufferedImage openBuffered(Path path) throws IOException {
        ImagePlus image = new Opener().openImage(path.toString());
        if (image == null) {
            throw new IOException("Could not open image: " + path);
        }
        int width = image.getWidth();
        int height = image.getHeight();
        BufferedImage buffered = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        ij.process.ImageProcessor processor = image.getProcessor().convertToColorProcessor();
        int[] rgb = new int[3];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                processor.getPixel(x, y, rgb);
                int packed = (rgb[0] << 16) | (rgb[1] << 8) | rgb[2];
                buffered.setRGB(x, y, packed);
            }
        }
        return buffered;
    }

    private static AnalysisFrame prepareAnalysisFrame(Path imagePath, BufferedImage source, Path outputDir) throws IOException {
        if (!isScanLike(source)) {
            throw new IllegalArgumentException(
                "This repository supports TIFF scans only. Photos and low-resolution images are intentionally not supported: " + imagePath
            );
        }

        BufferedImage analysis = source;
        Path analysisPath = imagePath;
        String sourceName = "full_image";
        Box analysisRoi = null;
        Box canvas = new Box(0, 0, source.getWidth(), source.getHeight());
        double scale = 1.0;

        if (isScanLike(source)) {
            int left = (int) Math.round(source.getWidth() * CONFIG.scanRoiLeftFraction());
            int top = (int) Math.round(source.getHeight() * CONFIG.scanRoiTopFraction());
            analysisRoi = new Box(left, top, source.getWidth(), source.getHeight());
            BufferedImage cropped = source.getSubimage(left, top, analysisRoi.width(), analysisRoi.height());

            int padX = Math.max(180, (int) Math.round(cropped.getWidth() * 0.16));
            int padY = Math.max(180, (int) Math.round(cropped.getHeight() * 0.10));
            BufferedImage padded = new BufferedImage(cropped.getWidth() + padX * 2, cropped.getHeight() + padY * 2, BufferedImage.TYPE_INT_RGB);
            int fill = medianRgb(cropped, Math.min(200, cropped.getWidth()), Math.min(200, cropped.getHeight()));
            for (int y = 0; y < padded.getHeight(); y++) {
                for (int x = 0; x < padded.getWidth(); x++) {
                    padded.setRGB(x, y, fill);
                }
            }
            padded.getGraphics().drawImage(cropped, padX, padY, null);
            canvas = new Box(left - padX, top - padY, left - padX + padded.getWidth(), top - padY + padded.getHeight());
            analysis = padded;
            sourceName = "scan_roi";

            int maxDim = Math.max(analysis.getWidth(), analysis.getHeight());
            if (maxDim > CONFIG.scanMaxDim()) {
                scale = CONFIG.scanMaxDim() / maxDim;
                int newW = Math.max(1, (int) Math.round(analysis.getWidth() * scale));
                int newH = Math.max(1, (int) Math.round(analysis.getHeight() * scale));
                BufferedImage resized = new BufferedImage(newW, newH, BufferedImage.TYPE_INT_RGB);
                resized.getGraphics().drawImage(analysis, 0, 0, newW, newH, null);
                analysis = resized;
                sourceName = "scan_roi_downscaled";
            }

            Path cropPath = outputDir.resolve("scan_plate_crop.png");
            ImageIO.write(analysis, "PNG", cropPath.toFile());
            saveScanRoiPreview(source, analysisRoi, outputDir.resolve("scan_plate_roi_preview.jpg"));
            analysisPath = cropPath;
        }

        return new AnalysisFrame(analysis, analysisPath, sourceName, analysisRoi, canvas, scale);
    }

    private static boolean isScanLike(BufferedImage image) {
        return image.getWidth() >= 5000 && image.getHeight() >= 5000;
    }

    private static void saveScanRoiPreview(BufferedImage source, Box roi, Path outputPath) throws IOException {
        int maxDim = 1600;
        double scale = Math.min(1.0, maxDim / (double) Math.max(source.getWidth(), source.getHeight()));
        int newW = (int) Math.round(source.getWidth() * scale);
        int newH = (int) Math.round(source.getHeight() * scale);
        BufferedImage preview = new BufferedImage(newW, newH, BufferedImage.TYPE_INT_RGB);
        preview.getGraphics().drawImage(source, 0, 0, newW, newH, null);
        java.awt.Graphics2D graphics = preview.createGraphics();
        graphics.setColor(new Color(255, 60, 60));
        graphics.drawRect(
            (int) Math.round(roi.left() * scale),
            (int) Math.round(roi.top() * scale),
            Math.max(1, (int) Math.round(roi.width() * scale)),
            Math.max(1, (int) Math.round(roi.height() * scale))
        );
        graphics.dispose();
        ImageIO.write(preview, "JPEG", outputPath.toFile());
    }

    private static int medianRgb(BufferedImage image, int sampleW, int sampleH) {
        int count = sampleW * sampleH;
        int[] rs = new int[count];
        int[] gs = new int[count];
        int[] bs = new int[count];
        int index = 0;
        for (int y = 0; y < sampleH; y++) {
            for (int x = 0; x < sampleW; x++) {
                int rgb = image.getRGB(x, y);
                rs[index] = (rgb >> 16) & 0xFF;
                gs[index] = (rgb >> 8) & 0xFF;
                bs[index] = rgb & 0xFF;
                index++;
            }
        }
        Arrays.sort(rs);
        Arrays.sort(gs);
        Arrays.sort(bs);
        int mid = count / 2;
        return (rs[mid] << 16) | (gs[mid] << 8) | bs[mid];
    }

    private static double[][] toGray(BufferedImage image) {
        double[][] gray = new double[image.getHeight()][image.getWidth()];
        for (int y = 0; y < image.getHeight(); y++) {
            for (int x = 0; x < image.getWidth(); x++) {
                int rgb = image.getRGB(x, y);
                int r = (rgb >> 16) & 0xFF;
                int g = (rgb >> 8) & 0xFF;
                int b = rgb & 0xFF;
                gray[y][x] = (0.2126 * r + 0.7152 * g + 0.0722 * b) / 255.0;
            }
        }
        return gray;
    }

    private static double[][] toPurple(BufferedImage image) {
        double[][] purple = new double[image.getHeight()][image.getWidth()];
        for (int y = 0; y < image.getHeight(); y++) {
            for (int x = 0; x < image.getWidth(); x++) {
                int rgb = image.getRGB(x, y);
                double r = ((rgb >> 16) & 0xFF) / 255.0;
                double g = ((rgb >> 8) & 0xFF) / 255.0;
                double b = (rgb & 0xFF) / 255.0;
                purple[y][x] = (r + b) * 0.5 - g;
            }
        }
        return purple;
    }

    private static Box detectPlateBox(double[][] gray) {
        double[][] smooth = gaussian(gray, 18.0);
        double[][] score = new double[gray.length][gray[0].length];
        double[] colProfile = new double[gray[0].length];
        double[] rowProfile = new double[gray.length];
        for (int y = 0; y < gray.length; y++) {
            for (int x = 0; x < gray[0].length; x++) {
                double value = Math.max(0.0, smooth[y][x] - gray[y][x]);
                score[y][x] = value;
                rowProfile[y] += value;
                colProfile[x] += value;
            }
        }
        colProfile = smoothProfile(colProfile, Math.max(11, gray[0].length / 40));
        rowProfile = smoothProfile(rowProfile, Math.max(11, gray.length / 40));

        int left = profileStart(colProfile);
        int right = profileEnd(colProfile);
        int top = profileStart(rowProfile);
        int bottom = profileEnd(rowProfile);
        int padX = Math.max(20, (right - left) / 16);
        int padY = Math.max(20, (bottom - top) / 16);
        left = Math.max(0, left - padX);
        top = Math.max(0, top - padY);
        right = Math.min(gray[0].length, right + padX);
        bottom = Math.min(gray.length, bottom + padY);
        return new Box(left, top, right, bottom);
    }

    private static double[] smoothProfile(double[] values, int window) {
        double[] out = new double[values.length];
        int radius = Math.max(1, window / 2);
        for (int i = 0; i < values.length; i++) {
            int start = Math.max(0, i - radius);
            int end = Math.min(values.length, i + radius + 1);
            double sum = 0.0;
            for (int j = start; j < end; j++) {
                sum += values[j];
            }
            out[i] = sum / (end - start);
        }
        return out;
    }

    private static int profileStart(double[] values) {
        double min = Arrays.stream(values).min().orElse(0.0);
        double max = Arrays.stream(values).max().orElse(0.0);
        double threshold = min + (max - min) * 0.34;
        for (int i = 0; i < values.length; i++) {
            if (values[i] >= threshold) {
                return i;
            }
        }
        return 0;
    }

    private static int profileEnd(double[] values) {
        double min = Arrays.stream(values).min().orElse(0.0);
        double max = Arrays.stream(values).max().orElse(0.0);
        double threshold = min + (max - min) * 0.34;
        for (int i = values.length - 1; i >= 0; i--) {
            if (values[i] >= threshold) {
                return i;
            }
        }
        return values.length;
    }

    private static List<Circle> createSeedCircles(Box plateBox, int width, int height) {
        boolean portrait = plateBox.height() > plateBox.width() * 1.1;
        List<Circle> seeds = new ArrayList<>();
        if (portrait) {
            double[] xs = new double[] {0.31, 0.69};
            double[] ys = new double[] {0.20, 0.50, 0.80};
            double radius = Math.min(plateBox.width() / 4.8, plateBox.height() / 7.0);
            int wellIndex = 1;
            for (double xf : xs) {
                for (double yf : ys) {
                    seeds.add(new Circle(plateBox.left() + plateBox.width() * xf, plateBox.top() + plateBox.height() * yf, radius, 1.0));
                    wellIndex++;
                }
            }
        } else {
            double[] xs = new double[] {0.19, 0.50, 0.81};
            double[] ys = new double[] {0.31, 0.69};
            double radius = Math.min(plateBox.width() / 7.0, plateBox.height() / 4.8);
            for (double yf : ys) {
                for (double xf : xs) {
                    seeds.add(new Circle(plateBox.left() + plateBox.width() * xf, plateBox.top() + plateBox.height() * yf, radius, 1.0));
                }
            }
        }
        return seeds;
    }

    private static List<Circle> createScanSeedCircles(int width, int height) {
        if (!CUSTOM_LAYOUT.isEmpty()) {
            List<Circle> seeds = new ArrayList<>();
            double minDim = Math.min(width, height);
            for (LayoutSeed seed : CUSTOM_LAYOUT) {
                seeds.add(
                    new Circle(
                        width * seed.xFraction(),
                        height * seed.yFraction(),
                        minDim * seed.radiusFraction(),
                        1.0
                    )
                );
            }
            return seeds;
        }

        boolean portrait = height > width * 1.1;
        List<Circle> seeds = new ArrayList<>();
        if (portrait) {
            double[] xs = new double[] {0.448, 0.732};
            double[] ys = new double[] {0.540, 0.688, 0.836};
            double radius = Math.min(width * 0.13, height * 0.072);
            int wellIndex = 1;
            for (double xf : xs) {
                for (double yf : ys) {
                    seeds.add(new Circle(width * xf, height * yf, radius, 1.0));
                    wellIndex++;
                }
            }
            return seeds;
        }

        double[] xs = new double[] {0.20, 0.50, 0.80};
        double[] ys = new double[] {0.31, 0.69};
        double radius = Math.min(width * 0.12, height * 0.16);
        for (double yf : ys) {
            for (double xf : xs) {
                seeds.add(new Circle(width * xf, height * yf, radius, 1.0));
            }
        }
        return seeds;
    }

    private static List<Circle> refineWells(Path analysisImagePath, Path seedCsv, Path outputCsv) throws Exception {
        WellMaskRefiner.main(new String[] {analysisImagePath.toString(), seedCsv.toString(), outputCsv.toString()});
        List<Circle> circles = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(outputCsv, StandardCharsets.UTF_8)) {
            reader.readLine();
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isBlank()) {
                    continue;
                }
                String[] parts = line.split(",");
                circles.add(
                    new Circle(
                        Double.parseDouble(parts[1]),
                        Double.parseDouble(parts[2]),
                        Double.parseDouble(parts[4]),
                        1.0
                    )
                );
            }
        }
        return circles;
    }

    private static void writeSeedCsv(Path path, List<Circle> seeds) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            writer.write("well_index,x,y,radius\n");
            for (int i = 0; i < seeds.size(); i++) {
                Circle seed = seeds.get(i);
                writer.write(String.format(LOCALE, "%d,%.4f,%.4f,%.4f%n", i + 1, seed.x(), seed.y(), seed.radius()));
            }
        }
    }

    private static List<WellRegion> buildWellRegions(List<Circle> circles, int width, int height) {
        List<WellRegion> regions = new ArrayList<>();
        for (int i = 0; i < circles.size(); i++) {
            Circle circle = circles.get(i);
            Circle rough = circle;
            int cropRadius = (int) Math.round(circle.radius() * CONFIG.wellCropScale());
            int left = Math.max(0, (int) Math.floor(circle.x() - cropRadius));
            int top = Math.max(0, (int) Math.floor(circle.y() - cropRadius));
            int right = Math.min(width, (int) Math.ceil(circle.x() + cropRadius));
            int bottom = Math.min(height, (int) Math.ceil(circle.y() + cropRadius));
            boolean[][] maskLocal = makeDiskMask(right - left, bottom - top, circle.x() - left, circle.y() - top, circle.radius());
            boolean[][] coreMask = erodeDisk(maskLocal, CONFIG.countMaskErosionPx());
            regions.add(new WellRegion(i + 1, rough, circle, new Box(left, top, right, bottom), maskLocal, coreMask));
        }
        return regions;
    }

    private static double[][] computeScoreMap(double[][] gray, double[][] purple, double meanRadius) {
        double sigma = Math.max(6.0, meanRadius * CONFIG.backgroundSigmaScale());
        double[][] grayBlur = gaussian(gray, sigma);
        double[][] purpleBlur = gaussian(purple, sigma);
        double[][] score = new double[gray.length][gray[0].length];
        for (int y = 0; y < gray.length; y++) {
            for (int x = 0; x < gray[0].length; x++) {
                double darkScore = grayBlur[y][x] - gray[y][x];
                double purpleScore = Math.max(0.0, purple[y][x] - purpleBlur[y][x]);
                score[y][x] = 0.7 * darkScore + 1.3 * purpleScore;
            }
        }
        return score;
    }

    private static WellOutput analyzeWell(double[][] scoreMap, double[][] gray, double[][] purple, WellRegion region, Metadata metadata, boolean scanLike) {
        double[][] scoreCrop = crop(scoreMap, region.cropBox());
        double[][] grayCrop = crop(gray, region.cropBox());
        double[][] purpleCrop = crop(purple, region.cropBox());
        boolean[][] mask = region.maskLocal();
        boolean[][] coreMask = region.coreMaskLocal();

        double[] insideValues = collectValues(scoreCrop, coreMask);
        double otsu = otsuThreshold(insideValues);
        double quantile = quantile(insideValues, CONFIG.scoreQuantile());
        double strict = Math.max(otsu, quantile);
        double loose = Math.min(otsu, quantile);
        double threshold = Math.max(loose, strict * CONFIG.thresholdRelaxation());
        double strongThreshold = scanLike ? Math.max(loose, threshold * CONFIG.scanThresholdFactor()) : threshold;

        boolean[][] strongBinary = new boolean[mask.length][mask[0].length];
        for (int y = 0; y < mask.length; y++) {
            for (int x = 0; x < mask[0].length; x++) {
                strongBinary[y][x] = mask[y][x] && scoreCrop[y][x] > strongThreshold;
            }
        }
        boolean[][] binary = strongBinary;
        binary = close(binary, 1);
        if (!scanLike) {
            binary = open(binary, 1);
        }
        binary = fillHoles(binary);
        binary = removeSmall(binary, CONFIG.minColonyArea());

        List<ConnectedComponent> mergedComponents = connectedComponents(binary);
        double[] mergedAreas = mergedComponents.stream().mapToDouble(component -> component.pixels().size()).toArray();
        double preSplitSingleArea = estimateSingleArea(mergedAreas, CONFIG.minColonyArea());
        boolean[][] splitMask = scanLike ? open(binary, 1) : binary;
        List<ConnectedComponent> resolvedComponents = selectiveSplitComponents(binary, splitMask, mergedComponents, preSplitSingleArea);

        List<RawParticle> rawParticles = analyzeParticles(resolvedComponents, scoreCrop, grayCrop, purpleCrop, mask, coreMask);
        double[] referenceAreas = rawParticles.stream()
            .filter(particle -> particle.edgeDistance() > 18.0 && particle.circularity() > 0.55)
            .mapToDouble(RawParticle::area)
            .toArray();
        double[] rawAreas = rawParticles.stream().mapToDouble(RawParticle::area).toArray();
        double singleArea = estimateSingleArea(referenceAreas.length > 0 ? referenceAreas : rawAreas, preSplitSingleArea);
        double[][] edgeDistance = distanceTransform(mask);

        List<SegmentRecord> kept = new ArrayList<>();
        WellAssignment assignment = resolveWellAssignment(metadata, region.wellIndex());
        String rowName = assignment.rowName();
        Double concentration = assignment.concentration();
        int replicate = assignment.replicate();
        int labelIndex = 0;
        double minAreaScale = scanLike ? CONFIG.minAreaScaleScan() : CONFIG.minAreaScalePhoto();

        for (RawParticle particle : rawParticles) {
            double localArea = localReferenceArea(particle, rawParticles, singleArea, region.refinedCircle().radius() * CONFIG.neighborRadiusScale());
            double minAcceptedArea = Math.max(CONFIG.minColonyArea(), localArea * minAreaScale);
            int support = supportNeighbors(particle, rawParticles, region.refinedCircle().radius() * CONFIG.neighborRadiusScale(), Math.max(CONFIG.minColonyArea() * 1.6, localArea * CONFIG.isolatedSupportAreaScale()));
            if (isRimArtifact(particle, localArea)) {
                continue;
            }
            if (scanLike
                && particle.area() < Math.max(CONFIG.minColonyArea() * TUNE_SCAN_SMALL_WEAK_AREA_MULTIPLIER, localArea * TUNE_SCAN_SMALL_WEAK_LOCAL_SCALE)
                && support <= TUNE_SCAN_SMALL_WEAK_MAX_SUPPORT
                && particle.meanGray() > TUNE_SCAN_SMALL_WEAK_MAX_GRAY
                && particle.meanPurple() < TUNE_SCAN_SMALL_WEAK_MIN_PURPLE
                && particle.meanScore() < Math.max(threshold * TUNE_SCAN_SMALL_WEAK_THRESHOLD_SCALE, TUNE_SCAN_SMALL_WEAK_MIN_SCORE)) {
                continue;
            }
            if (particle.area() < minAcceptedArea) {
                boolean borderlineStrong =
                    particle.area() >= minAcceptedArea * TUNE_BORDERLINE_AREA_SCALE
                    && particle.solidity() > TUNE_BORDERLINE_MIN_SOLIDITY
                    && particle.circularity() > TUNE_BORDERLINE_MIN_CIRCULARITY
                    && (
                        particle.meanGray() < TUNE_BORDERLINE_MAX_GRAY
                        || particle.meanPurple() > TUNE_BORDERLINE_MIN_PURPLE
                        || particle.meanScore() > Math.max(threshold * TUNE_BORDERLINE_THRESHOLD_SCALE, TUNE_BORDERLINE_MIN_SCORE)
                    );
                if (!borderlineStrong) {
                    continue;
                }
            }
            if (particle.area() < Math.max(CONFIG.minColonyArea() * 1.5, localArea * CONFIG.isolatedSmallAreaScale())
                && support == 0
                && particle.meanGray() > 0.11
                && particle.meanPurple() < 0.12) {
                continue;
            }
            kept.add(toSegmentRecord(particle, region, rowName, concentration, replicate, assignment.extraFields(), labelIndex++));
        }

        boolean[][] countedLocal = new boolean[mask.length][mask[0].length];
        for (SegmentRecord segment : kept) {
            paintRoiMask(countedLocal, segment.roiGlobal(), region.cropBox().left(), region.cropBox().top());
        }

        if (scanLike) {
            double weakThreshold = Math.max(loose, threshold * TUNE_SCAN_RESCUE_THRESHOLD_FACTOR);
            if (weakThreshold < strongThreshold) {
                boolean[][] weakBinary = new boolean[mask.length][mask[0].length];
                for (int y = 0; y < mask.length; y++) {
                    for (int x = 0; x < mask[0].length; x++) {
                        weakBinary[y][x] = mask[y][x] && scoreCrop[y][x] > weakThreshold;
                    }
                }
                weakBinary = close(weakBinary, 1);
                weakBinary = fillHoles(weakBinary);
                weakBinary = removeSmall(weakBinary, Math.max(6, CONFIG.minColonyArea() / 2));

                List<ConnectedComponent> weakMerged = connectedComponents(weakBinary);
                boolean[][] weakSplitMask = open(weakBinary, 1);
                List<ConnectedComponent> weakResolved = selectiveSplitComponents(weakBinary, weakSplitMask, weakMerged, singleArea);
                List<RawParticle> weakParticles = analyzeParticles(weakResolved, scoreCrop, grayCrop, purpleCrop, mask, coreMask);

                for (RawParticle particle : weakParticles) {
                    double overlapFraction = roiOverlapFraction(countedLocal, particle.roiLocal());
                    double novelFraction = 1.0 - overlapFraction;
                    if (overlapFraction > TUNE_SCAN_RESCUE_MAX_OVERLAP || novelFraction < TUNE_SCAN_RESCUE_MIN_NEW_FRACTION) {
                        continue;
                    }

                    double localArea = localReferenceArea(particle, rawParticles, singleArea, region.refinedCircle().radius() * CONFIG.neighborRadiusScale());
                    double rescueMinArea = Math.max(CONFIG.minColonyArea() * 2.0, localArea * TUNE_SCAN_RESCUE_MIN_AREA_SCALE);
                    boolean rescueCandidate =
                        particle.area() >= rescueMinArea
                        && particle.solidity() >= TUNE_SCAN_RESCUE_MIN_SOLIDITY
                        && (
                            particle.circularity() >= TUNE_SCAN_RESCUE_MIN_CIRCULARITY
                            || particle.meanScore() >= TUNE_SCAN_RESCUE_MIN_SCORE
                        )
                        && (
                            particle.meanGray() <= TUNE_SCAN_RESCUE_MAX_GRAY
                            || particle.meanPurple() >= TUNE_SCAN_RESCUE_MIN_PURPLE
                            || particle.meanScore() >= TUNE_SCAN_RESCUE_MIN_SCORE
                        );
                    if (!rescueCandidate || isRimArtifact(particle, localArea)) {
                        continue;
                    }

                    SegmentRecord rescued = toSegmentRecord(particle, region, rowName, concentration, replicate, assignment.extraFields(), labelIndex++);
                    kept.add(rescued);
                    paintRoiMask(countedLocal, rescued.roiGlobal(), region.cropBox().left(), region.cropBox().top());
                }

                mergeCompactWeakBlobs:
                {
                    Set<Integer> removed = new HashSet<>();
                    List<SegmentRecord> mergedBack = new ArrayList<>();
                    List<RawParticle> mergeParticles = new ArrayList<>(weakParticles);
                    mergeParticles.sort(Comparator.comparingInt(RawParticle::area).reversed());
                    for (RawParticle particle : mergeParticles) {
                        if (particle.solidity() < TUNE_SCAN_MERGE_MIN_SOLIDITY
                            || particle.circularity() < TUNE_SCAN_MERGE_MIN_CIRCULARITY
                            || particle.bboxAspectRatio() > TUNE_SCAN_MERGE_MAX_ASPECT) {
                            continue;
                        }

                        List<Integer> members = new ArrayList<>();
                        int coveredArea = 0;
                        for (int index = 0; index < kept.size(); index++) {
                            if (removed.contains(index)) {
                                continue;
                            }
                            SegmentRecord segment = kept.get(index);
                            int localX = (int) Math.round(segment.centroidX() - region.cropBox().left());
                            int localY = (int) Math.round(segment.centroidY() - region.cropBox().top());
                            if (roiContainsPoint(particle.roiLocal(), localX, localY)) {
                                members.add(index);
                                coveredArea += segment.areaPx();
                            }
                        }

                        if (members.size() < TUNE_SCAN_MERGE_MIN_SEGMENTS
                            || coveredArea < particle.area() * TUNE_SCAN_MERGE_MIN_COVERAGE) {
                            continue;
                        }

                        for (int index : members) {
                            removed.add(index);
                        }
                        mergedBack.add(toSegmentRecord(particle, region, rowName, concentration, replicate, assignment.extraFields(), labelIndex++));
                    }

                    if (!removed.isEmpty()) {
                        List<SegmentRecord> merged = new ArrayList<>();
                        for (int index = 0; index < kept.size(); index++) {
                            if (!removed.contains(index)) {
                                merged.add(kept.get(index));
                            }
                        }
                        merged.addAll(mergedBack);
                        kept = merged;

                        countedLocal = new boolean[mask.length][mask[0].length];
                        for (SegmentRecord segment : kept) {
                            paintRoiMask(countedLocal, segment.roiGlobal(), region.cropBox().left(), region.cropBox().top());
                        }
                    }
                }
            }
        }

        kept = reindexSegments(kept);

        WellSummary summary = new WellSummary(
            region.wellIndex(),
            rowName,
            replicate,
            concentration,
            assignment.extraFields(),
            rawParticles.size(),
            kept.size(),
            singleArea,
            threshold,
            countTrue(mask),
            countTrue(coreMask)
        );
        return new WellOutput(kept, summary, countedLocal);
    }

    private static List<RawParticle> analyzeParticles(
        List<ConnectedComponent> components,
        double[][] score,
        double[][] gray,
        double[][] purple,
        boolean[][] innerMask,
        boolean[][] coreMask
    ) {
        double[][] edgeDistance = distanceTransform(innerMask);
        List<RawParticle> particles = new ArrayList<>();
        for (int i = 0; i < components.size(); i++) {
            ConnectedComponent component = components.get(i);
            Roi roi = roiFromComponent(component);
            ParticleStats stats = particleStats(roi, score, gray, purple, edgeDistance, coreMask);
            double area = component.pixels().size();
            double perimeter = stats.perimeter;
            double circularity = stats.circularity;
            double solidity = stats.solidity;
            double aspectRatio = stats.aspectRatio;
            double eccentricity = aspectRatio <= 1.0 ? 0.0 : Math.sqrt(1.0 - 1.0 / (aspectRatio * aspectRatio));
            particles.add(
                new RawParticle(
                    i,
                    roi,
                    (int) Math.round(area),
                    perimeter,
                    circularity,
                    solidity,
                    eccentricity,
                    stats.centroidX,
                    stats.centroidY,
                    stats.meanScore,
                    stats.meanGray,
                    stats.meanPurple,
                    stats.edgeDistance,
                    stats.coreOverlap,
                    roi.getBounds()
                )
            );
        }
        return particles;
    }

    private static List<ConnectedComponent> selectiveSplitComponents(
        boolean[][] parentMask,
        boolean[][] splitMask,
        List<ConnectedComponent> mergedComponents,
        double referenceArea
    ) {
        if (mergedComponents.isEmpty()) {
            return mergedComponents;
        }

        ByteProcessor watershedInput = toByteProcessor(splitMask);
        new EDM().toWatershed(watershedInput);
        boolean[][] watershedMask = fromByteProcessor(watershedInput);
        watershedMask = removeSmall(watershedMask, CONFIG.minColonyArea());
        List<ConnectedComponent> splitComponents = connectedComponents(watershedMask);
        if (splitComponents.isEmpty()) {
            return mergedComponents;
        }

        int[][] parentLabels = componentLabelMap(mergedComponents, parentMask.length, parentMask[0].length);
        Map<Integer, List<ConnectedComponent>> splitByParent = new HashMap<>();
        for (ConnectedComponent component : splitComponents) {
            int[] first = component.pixels().getFirst();
            int parentId = parentLabels[first[0]][first[1]];
            if (parentId <= 0) {
                continue;
            }
            splitByParent.computeIfAbsent(parentId, ignored -> new ArrayList<>()).add(component);
        }

        double splitTriggerArea = Math.max(referenceArea * TUNE_SPLIT_TRIGGER_AREA_SCALE, CONFIG.minColonyArea() * 4.0);
        double minSplitArea = Math.max(referenceArea * TUNE_SPLIT_MIN_CHILD_AREA_SCALE, CONFIG.minColonyArea() * 2.0);
        List<ConnectedComponent> selected = new ArrayList<>();

        for (int i = 0; i < mergedComponents.size(); i++) {
            ConnectedComponent parent = mergedComponents.get(i);
            double parentArea = parent.pixels().size();
            List<ConnectedComponent> children = splitByParent.getOrDefault(i + 1, List.of());
            if (children.size() < 2 || parentArea < splitTriggerArea) {
                selected.add(parent);
                continue;
            }

            ComponentShape parentShape = componentShape(parent);
            boolean compactParent =
                parentShape.solidity() >= TUNE_SPLIT_COMPACT_MIN_SOLIDITY
                && parentShape.circularity() >= TUNE_SPLIT_COMPACT_MIN_CIRCULARITY
                && parentShape.aspectRatio() <= TUNE_SPLIT_COMPACT_MAX_ASPECT;
            if (compactParent && parentArea < referenceArea * TUNE_SPLIT_COMPACT_TRIGGER_AREA_SCALE) {
                selected.add(parent);
                continue;
            }

            List<ConnectedComponent> validChildren = new ArrayList<>();
            int validCoverage = 0;
            children.sort(Comparator.comparingInt((ConnectedComponent component) -> component.pixels().size()).reversed());
            for (ConnectedComponent child : children) {
                if (child.pixels().size() >= minSplitArea) {
                    validChildren.add(child);
                    validCoverage += child.pixels().size();
                }
            }

            double childSpread = maxCentroidDistance(validChildren);
            double equivalentDiameter = 2.0 * Math.sqrt(Math.max(referenceArea, 1.0) / Math.PI);
            boolean wellSeparatedChildren = childSpread >= equivalentDiameter * TUNE_SPLIT_MIN_CENTROID_DISTANCE_SCALE;

            if (validChildren.size() >= 2
                && validCoverage >= parentArea * TUNE_SPLIT_VALID_COVERAGE
                && (!compactParent || wellSeparatedChildren || validChildren.size() >= 3)) {
                selected.addAll(validChildren);
            } else {
                selected.add(parent);
            }
        }
        return selected;
    }

    private static double maxCentroidDistance(List<ConnectedComponent> components) {
        double maxDistance = 0.0;
        for (int i = 0; i < components.size(); i++) {
            ComponentShape a = componentShape(components.get(i));
            for (int j = i + 1; j < components.size(); j++) {
                ComponentShape b = componentShape(components.get(j));
                maxDistance = Math.max(maxDistance, Math.hypot(a.centroidX() - b.centroidX(), a.centroidY() - b.centroidY()));
            }
        }
        return maxDistance;
    }

    private static ComponentShape componentShape(ConnectedComponent component) {
        Roi roi = roiFromComponent(component);
        double area = component.pixels().size();
        double perimeter = roi.getLength();
        double circularity = perimeter <= 0.0 ? 0.0 : (4.0 * Math.PI * area / (perimeter * perimeter));
        List<double[]> points = new ArrayList<>(component.pixels().size());
        double sumX = 0.0;
        double sumY = 0.0;
        for (int[] point : component.pixels()) {
            points.add(new double[] {point[1], point[0]});
            sumX += point[1];
            sumY += point[0];
        }
        double solidity = convexHullArea(points) <= 0.0 ? 1.0 : Math.min(1.0, area / convexHullArea(points));
        Rectangle bounds = component.bounds();
        double width = Math.max(1.0, bounds.width);
        double height = Math.max(1.0, bounds.height);
        double aspectRatio = Math.max(width, height) / Math.min(width, height);
        return new ComponentShape(
            circularity,
            solidity,
            aspectRatio,
            sumX / Math.max(area, 1.0),
            sumY / Math.max(area, 1.0)
        );
    }

    private static int[][] componentLabelMap(List<ConnectedComponent> components, int height, int width) {
        int[][] labels = new int[height][width];
        for (int i = 0; i < components.size(); i++) {
            for (int[] point : components.get(i).pixels()) {
                labels[point[0]][point[1]] = i + 1;
            }
        }
        return labels;
    }

    private static Roi roiFromComponent(ConnectedComponent component) {
        Rectangle bounds = component.bounds();
        ByteProcessor mask = new ByteProcessor(bounds.width, bounds.height);
        for (int[] point : component.pixels()) {
            mask.set(point[1] - bounds.x, point[0] - bounds.y, 255);
        }

        int[] seed = component.pixels().getFirst();
        Wand wand = new Wand(mask);
        wand.autoOutline(seed[1] - bounds.x, seed[0] - bounds.y, 255, 255, Wand.EIGHT_CONNECTED);
        if (wand.npoints <= 0) {
            return new Roi(bounds.x, bounds.y, bounds.width, bounds.height);
        }

        int[] xs = new int[wand.npoints];
        int[] ys = new int[wand.npoints];
        for (int i = 0; i < wand.npoints; i++) {
            xs[i] = wand.xpoints[i] + bounds.x;
            ys[i] = wand.ypoints[i] + bounds.y;
        }
        return new PolygonRoi(xs, ys, wand.npoints, Roi.TRACED_ROI);
    }

    private static boolean isRimArtifact(RawParticle particle, double localAreaRef) {
        boolean isBoundary = particle.edgeDistance() < 20.0 || particle.coreOverlap() < 0.45;
        boolean touchesRim = particle.edgeDistance() < 4.0;
        boolean strongCandidate =
            particle.meanGray() < 0.12
            || (particle.meanPurple() > 0.12 && particle.meanScore() > 0.08)
            || (particle.meanScore() > 0.11 && particle.solidity() > 0.80);

        boolean reject = particle.edgeDistance() <= CONFIG.outerRejectBandPx() && !strongCandidate;
        if (isBoundary) {
            if (touchesRim && particle.meanGray() > 0.24 && particle.meanPurple() < 0.09) {
                reject = true;
            } else if (touchesRim && particle.meanGray() > 0.19 && particle.meanPurple() < 0.055) {
                reject = true;
            } else if (touchesRim && particle.eccentricity() > 0.95 && particle.circularity() < 0.55) {
                reject = true;
            } else if (particle.edgeDistance() < 26.0 && particle.meanGray() > 0.34) {
                reject = true;
            } else if (particle.edgeDistance() < 28.0 && particle.meanGray() > 0.25 && particle.meanPurple() < 0.045) {
                reject = true;
            }
            if (particle.coreOverlap() < 0.12 && particle.meanScore() < 0.08) {
                reject = true;
            } else if (particle.coreOverlap() < 0.18 && particle.meanPurple() < 0.065) {
                reject = true;
            } else if (particle.coreOverlap() < 0.25 && particle.meanGray() > 0.17 && particle.meanPurple() < 0.085) {
                reject = true;
            } else if (particle.edgeDistance() < 16.0 && particle.meanGray() > 0.22 && particle.meanPurple() < 0.08) {
                reject = true;
            } else if (particle.edgeDistance() < 18.0 && particle.meanGray() > 0.24 && particle.meanPurple() < 0.085) {
                reject = true;
            } else if (particle.edgeDistance() < 20.0 && particle.meanGray() > 0.32 && particle.meanPurple() < 0.10) {
                reject = true;
            }
            if (strongCandidate && particle.meanPurple() > 0.10 && particle.meanGray() < 0.15) {
                reject = false;
            } else if (strongCandidate && particle.meanScore() > 0.12 && particle.solidity() > 0.82) {
                reject = false;
            }
            if (particle.area() > Math.max(localAreaRef * 3.8, 260.0) && particle.circularity() < 0.48 && particle.coreOverlap() < 0.55) {
                reject = true;
            }
            if (particle.edgeDistance() < 10.0
                && particle.bboxAspectRatio() > 1.8
                && particle.circularity() < 0.78
                && particle.meanGray() > 0.18) {
                reject = true;
            } else if (particle.edgeDistance() < 15.0
                && particle.bboxAspectRatio() >= 1.9
                && particle.circularity() < 0.60
                && particle.meanGray() > 0.20) {
                reject = true;
            } else if (particle.edgeDistance() < 18.0
                && particle.bboxAspectRatio() > 2.2
                && particle.circularity() < 0.72) {
                reject = true;
            } else if (particle.edgeDistance() < 12.0
                && particle.bboxAspectRatio() > 1.6
                && particle.bboxExtent() < 0.68
                && particle.meanPurple() < 0.055
                && particle.meanGray() > 0.20) {
                reject = true;
            }
            if (particle.edgeDistance() < 12.0
                && particle.area() < Math.max(localAreaRef * 0.90, 70.0)
                && (particle.circularity() < 0.88 || particle.eccentricity() > 0.80)
                && particle.meanGray() > 0.16) {
                reject = true;
            } else if (particle.edgeDistance() < 18.0
                && particle.area() < Math.max(localAreaRef * 0.75, 60.0)
                && particle.circularity() < 0.72
                && particle.eccentricity() > 0.90) {
                reject = true;
            }
        }
        return reject;
    }

    private static double localReferenceArea(RawParticle target, List<RawParticle> particles, double fallback, double radius) {
        List<Double> nearby = new ArrayList<>();
        for (RawParticle particle : particles) {
            if (particle.id() == target.id()) {
                continue;
            }
            if (Math.hypot(target.centroidX() - particle.centroidX(), target.centroidY() - particle.centroidY()) <= radius) {
                nearby.add((double) particle.area());
            }
        }
        if (nearby.size() >= 4) {
            return estimateSingleArea(nearby.stream().mapToDouble(Double::doubleValue).toArray(), fallback);
        }
        return fallback;
    }

    private static int supportNeighbors(RawParticle target, List<RawParticle> particles, double radius, double minArea) {
        int count = 0;
        for (RawParticle particle : particles) {
            if (particle.id() == target.id() || particle.area() < minArea) {
                continue;
            }
            if (Math.hypot(target.centroidX() - particle.centroidX(), target.centroidY() - particle.centroidY()) <= radius) {
                count++;
            }
        }
        return count;
    }

    private static double estimateSingleArea(double[] values, double fallback) {
        if (values.length == 0) {
            return fallback;
        }
        Arrays.sort(values);
        if (values.length < 5) {
            return values[values.length / 2];
        }
        double low = quantile(values, TUNE_SINGLE_AREA_LOW_QUANTILE);
        double high = quantile(values, TUNE_SINGLE_AREA_HIGH_QUANTILE);
        List<Double> core = new ArrayList<>();
        for (double value : values) {
            if (value >= low && value <= high) {
                core.add(value);
            }
        }
        if (core.isEmpty()) {
            return values[values.length / 2];
        }
        double[] coreArray = core.stream().mapToDouble(Double::doubleValue).toArray();
        Arrays.sort(coreArray);
        return quantile(coreArray, TUNE_SINGLE_AREA_TARGET_QUANTILE);
    }

    private static void saveWellDebug(BufferedImage analysisImage, double[][] scoreMap, WellRegion region, WellOutput output, Path outputDir) throws IOException {
        BufferedImage crop = analysisImage.getSubimage(region.cropBox().left(), region.cropBox().top(), region.cropBox().width(), region.cropBox().height());
        ImagePlus imp = new ImagePlus("well", crop);
        Overlay overlay = new Overlay();

        OvalRoi well = new OvalRoi(
            region.refinedCircle().x() - region.refinedCircle().radius() - region.cropBox().left(),
            region.refinedCircle().y() - region.refinedCircle().radius() - region.cropBox().top(),
            region.refinedCircle().radius() * 2.0,
            region.refinedCircle().radius() * 2.0
        );
        well.setStrokeColor(new Color(255, 210, 0));
        overlay.add(well);

        for (SegmentRecord segment : output.segments()) {
            Roi shifted = shiftRoi(segment.roiGlobal(), -region.cropBox().left(), -region.cropBox().top());
            shifted.setStrokeColor(Color.RED);
            overlay.add(shifted);
            TextRoi text = new TextRoi((float) (segment.centroidX() - region.cropBox().left() + 4), (float) (segment.centroidY() - region.cropBox().top() - 4), Integer.toString(segment.labelIndex()), LABEL_FONT);
            text.setStrokeColor(Color.YELLOW);
            overlay.add(text);
        }
        TextRoi header = new TextRoi(12, 12, "W" + region.wellIndex() + ": " + output.summary().countedColonies(), WELL_FONT);
        header.setStrokeColor(new Color(255, 210, 0));
        overlay.add(header);

        imp.setOverlay(overlay);
        ImagePlus flat = imp.flatten();
        new FileSaver(flat).saveAsPng(outputDir.resolve(String.format(LOCALE, "well_%02d_overlay.png", region.wellIndex())).toString());

        ByteProcessor mask = toByteProcessor(output.countedMaskLocal());
        new FileSaver(new ImagePlus("mask", mask)).saveAsPng(outputDir.resolve(String.format(LOCALE, "well_%02d_mask.png", region.wellIndex())).toString());

        double[][] scoreCrop = crop(scoreMap, region.cropBox());
        FloatProcessor scoreProcessor = toNormalizedProcessor(scoreCrop);
        new FileSaver(new ImagePlus("score", scoreProcessor.convertToByteProcessor())).saveAsPng(outputDir.resolve(String.format(LOCALE, "well_%02d_score.png", region.wellIndex())).toString());
    }

    private static void savePlateOverlay(BufferedImage analysisImage, List<WellRegion> regions, List<SegmentRecord> segments, List<WellSummary> wellSummaries, Path outputPath) {
        ImagePlus imp = new ImagePlus("plate", analysisImage);
        Overlay overlay = new Overlay();

        for (SegmentRecord segment : segments) {
            Roi roi = (Roi) segment.roiGlobal().clone();
            roi.setStrokeColor(Color.RED);
            overlay.add(roi);
            TextRoi text = new TextRoi((float) (segment.centroidX() + 4), (float) (segment.centroidY() - 4), Integer.toString(segment.labelIndex()), LABEL_FONT);
            text.setStrokeColor(Color.YELLOW);
            overlay.add(text);
        }

        for (int i = 0; i < regions.size(); i++) {
            WellRegion region = regions.get(i);
            WellSummary summary = wellSummaries.get(i);
            OvalRoi well = new OvalRoi(
                region.refinedCircle().x() - region.refinedCircle().radius(),
                region.refinedCircle().y() - region.refinedCircle().radius(),
                region.refinedCircle().radius() * 2.0,
                region.refinedCircle().radius() * 2.0
            );
            well.setStrokeColor(new Color(255, 210, 0));
            overlay.add(well);
            TextRoi text = new TextRoi((float) (well.getBounds().x), (float) (well.getBounds().y - 18), "W" + summary.wellIndex() + ": " + summary.countedColonies(), WELL_FONT);
            text.setStrokeColor(new Color(255, 210, 0));
            overlay.add(text);
        }

        imp.setOverlay(overlay);
        new FileSaver(imp.flatten()).saveAsPng(outputPath.toString());
    }

    private static void saveColoniesOnly(BufferedImage analysisImage, boolean[][] globalMask, Path outputPath) throws IOException {
        BufferedImage rgba = new BufferedImage(analysisImage.getWidth(), analysisImage.getHeight(), BufferedImage.TYPE_INT_ARGB);
        for (int y = 0; y < analysisImage.getHeight(); y++) {
            for (int x = 0; x < analysisImage.getWidth(); x++) {
                int rgb = analysisImage.getRGB(x, y) & 0x00FFFFFF;
                int alpha = globalMask[y][x] ? 0xFF : 0x00;
                rgba.setRGB(x, y, (alpha << 24) | rgb);
            }
        }
        ImageIO.write(rgba, "PNG", outputPath.toFile());
    }

    private static void mergeLocalMask(boolean[][] globalMask, boolean[][] localMask, Box cropBox) {
        for (int y = 0; y < localMask.length; y++) {
            for (int x = 0; x < localMask[0].length; x++) {
                if (localMask[y][x]) {
                    globalMask[cropBox.top() + y][cropBox.left() + x] = true;
                }
            }
        }
    }

    private static SegmentRecord toSegmentRecord(
        RawParticle particle,
        WellRegion region,
        String rowName,
        Double concentration,
        int replicate,
        Map<String, Object> extraFields,
        int labelIndex
    ) {
        Roi shifted = shiftRoi(particle.roiLocal(), region.cropBox().left(), region.cropBox().top());
        Rectangle bbox = shifted.getBounds();
        return new SegmentRecord(
            region.wellIndex(),
            rowName,
            concentration,
            replicate,
            new LinkedHashMap<>(extraFields),
            particle.area(),
            particle.eccentricity(),
            particle.circularity(),
            particle.solidity(),
            particle.edgeDistance(),
            particle.coreOverlap(),
            particle.meanGray(),
            particle.meanPurple(),
            particle.centroidX() + region.cropBox().left(),
            particle.centroidY() + region.cropBox().top(),
            bbox,
            shifted,
            labelIndex
        );
    }

    private static List<SegmentRecord> reindexSegments(List<SegmentRecord> segments) {
        List<SegmentRecord> relabeled = new ArrayList<>(segments.size());
        for (int i = 0; i < segments.size(); i++) {
            SegmentRecord segment = segments.get(i);
            relabeled.add(
                new SegmentRecord(
                    segment.wellIndex(),
                    segment.rowName(),
                    segment.concentration(),
                    segment.replicate(),
                    new LinkedHashMap<>(segment.extraFields()),
                    segment.areaPx(),
                    segment.eccentricity(),
                    segment.circularity(),
                    segment.solidity(),
                    segment.edgeDistancePx(),
                    segment.coreOverlap(),
                    segment.meanGray(),
                    segment.meanPurple(),
                    segment.centroidX(),
                    segment.centroidY(),
                    segment.bbox(),
                    segment.roiGlobal(),
                    i
                )
            );
        }
        return relabeled;
    }

    private static void paintRoiMask(boolean[][] mask, Roi roiGlobal, int offsetX, int offsetY) {
        Rectangle bounds = roiGlobal.getBounds();
        ByteProcessor roiMask = roiMask(roiGlobal);
        for (int y = 0; y < bounds.height; y++) {
            for (int x = 0; x < bounds.width; x++) {
                if (roiMask.get(x, y) == 0) {
                    continue;
                }
                int yy = bounds.y - offsetY + y;
                int xx = bounds.x - offsetX + x;
                if (yy >= 0 && xx >= 0 && yy < mask.length && xx < mask[0].length) {
                    mask[yy][xx] = true;
                }
            }
        }
    }

    private static double roiOverlapFraction(boolean[][] mask, Roi roiLocal) {
        Rectangle bounds = roiLocal.getBounds();
        ByteProcessor roiMask = roiMask(roiLocal);
        int total = 0;
        int overlap = 0;
        for (int y = 0; y < bounds.height; y++) {
            for (int x = 0; x < bounds.width; x++) {
                if (roiMask.get(x, y) == 0) {
                    continue;
                }
                int yy = bounds.y + y;
                int xx = bounds.x + x;
                if (yy < 0 || xx < 0 || yy >= mask.length || xx >= mask[0].length) {
                    continue;
                }
                total++;
                if (mask[yy][xx]) {
                    overlap++;
                }
            }
        }
        return total == 0 ? 0.0 : overlap / (double) total;
    }

    private static boolean roiContainsPoint(Roi roiLocal, int x, int y) {
        Rectangle bounds = roiLocal.getBounds();
        if (x < bounds.x || y < bounds.y || x >= bounds.x + bounds.width || y >= bounds.y + bounds.height) {
            return false;
        }
        ByteProcessor roiMask = roiMask(roiLocal);
        return roiMask.get(x - bounds.x, y - bounds.y) != 0;
    }

    private static ByteProcessor roiMask(Roi roi) {
        Rectangle bounds = roi.getBounds();
        ByteProcessor mask = new ByteProcessor(bounds.width, bounds.height);
        ByteProcessor roiMask = (ByteProcessor) roi.getMask();
        if (roiMask == null) {
            mask.setValue(255);
            mask.fill();
            return mask;
        }
        mask.insert(roiMask, 0, 0);
        return mask;
    }

    private static Roi shiftRoi(Roi roi, int dx, int dy) {
        FloatPolygon polygon = roi.getFloatPolygon();
        if (polygon.npoints == 0) {
            Roi shifted = (Roi) roi.clone();
            shifted.setLocation(roi.getBounds().x + dx, roi.getBounds().y + dy);
            return shifted;
        }
        float[] xs = new float[polygon.npoints];
        float[] ys = new float[polygon.npoints];
        for (int i = 0; i < polygon.npoints; i++) {
            xs[i] = polygon.xpoints[i] + dx;
            ys[i] = polygon.ypoints[i] + dy;
        }
        return new PolygonRoi(new FloatPolygon(xs, ys), Roi.POLYGON);
    }

    private static ParticleStats particleStats(
        Roi roi,
        double[][] score,
        double[][] gray,
        double[][] purple,
        double[][] edgeDistance,
        boolean[][] coreMask
    ) {
        Rectangle bounds = roi.getBounds();
        ByteProcessor mask = roiMask(roi);
        double area = 0.0;
        double sumX = 0.0;
        double sumY = 0.0;
        double sumScore = 0.0;
        double sumGray = 0.0;
        double sumPurple = 0.0;
        double sumCore = 0.0;
        double minEdge = Double.POSITIVE_INFINITY;
        List<double[]> points = new ArrayList<>();
        for (int y = 0; y < bounds.height; y++) {
            for (int x = 0; x < bounds.width; x++) {
                if (mask.get(x, y) == 0) {
                    continue;
                }
                int xx = bounds.x + x;
                int yy = bounds.y + y;
                area += 1.0;
                sumX += xx;
                sumY += yy;
                sumScore += score[yy][xx];
                sumGray += gray[yy][xx];
                sumPurple += purple[yy][xx];
                sumCore += coreMask[yy][xx] ? 1.0 : 0.0;
                minEdge = Math.min(minEdge, edgeDistance[yy][xx]);
                points.add(new double[] {xx, yy});
            }
        }
        double centroidX = sumX / Math.max(area, 1.0);
        double centroidY = sumY / Math.max(area, 1.0);
        double mu20 = 0.0;
        double mu02 = 0.0;
        double mu11 = 0.0;
        for (double[] point : points) {
            double dx = point[0] - centroidX;
            double dy = point[1] - centroidY;
            mu20 += dx * dx;
            mu02 += dy * dy;
            mu11 += dx * dy;
        }
        double trace = mu20 + mu02;
        double det = mu20 * mu02 - mu11 * mu11;
        double term = Math.sqrt(Math.max(0.0, trace * trace * 0.25 - det));
        double lambda1 = Math.max(1e-9, trace * 0.5 + term);
        double lambda2 = Math.max(1e-9, trace * 0.5 - term);
        double aspectRatio = Math.sqrt(lambda1 / lambda2);
        double perimeter = roi.getLength();
        double circularity = perimeter <= 0.0 ? 0.0 : (4.0 * Math.PI * area / (perimeter * perimeter));
        double solidity = convexHullArea(points) <= 0.0 ? 1.0 : Math.min(1.0, area / convexHullArea(points));
        return new ParticleStats(
            centroidX,
            centroidY,
            sumScore / Math.max(area, 1.0),
            sumGray / Math.max(area, 1.0),
            sumPurple / Math.max(area, 1.0),
            minEdge,
            sumCore / Math.max(area, 1.0),
            perimeter,
            circularity,
            solidity,
            aspectRatio
        );
    }

    private record ParticleStats(
        double centroidX,
        double centroidY,
        double meanScore,
        double meanGray,
        double meanPurple,
        double edgeDistance,
        double coreOverlap,
        double perimeter,
        double circularity,
        double solidity,
        double aspectRatio
    ) {}

    private static double convexHullArea(List<double[]> points) {
        if (points.size() < 3) {
            return points.size();
        }
        points.sort((a, b) -> a[0] == b[0] ? Double.compare(a[1], b[1]) : Double.compare(a[0], b[0]));
        List<double[]> lower = new ArrayList<>();
        for (double[] point : points) {
            while (lower.size() >= 2 && cross(lower.get(lower.size() - 2), lower.get(lower.size() - 1), point) <= 0.0) {
                lower.remove(lower.size() - 1);
            }
            lower.add(point);
        }
        List<double[]> upper = new ArrayList<>();
        for (int i = points.size() - 1; i >= 0; i--) {
            double[] point = points.get(i);
            while (upper.size() >= 2 && cross(upper.get(upper.size() - 2), upper.get(upper.size() - 1), point) <= 0.0) {
                upper.remove(upper.size() - 1);
            }
            upper.add(point);
        }
        lower.remove(lower.size() - 1);
        upper.remove(upper.size() - 1);
        List<double[]> hull = new ArrayList<>(lower);
        hull.addAll(upper);
        double area = 0.0;
        for (int i = 0; i < hull.size(); i++) {
            double[] a = hull.get(i);
            double[] b = hull.get((i + 1) % hull.size());
            area += a[0] * b[1] - b[0] * a[1];
        }
        return Math.abs(area) * 0.5;
    }

    private static double cross(double[] a, double[] b, double[] c) {
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    }

    private static double safeTable(ResultsTable table, String column, int row, double fallback) {
        if (table.columnExists(column)) {
            return table.getValue(column, row);
        }
        return fallback;
    }

    private static FloatProcessor toNormalizedProcessor(double[][] values) {
        int width = values[0].length;
        int height = values.length;
        float[] pixels = new float[width * height];
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for (double[] row : values) {
            for (double value : row) {
                min = Math.min(min, value);
                max = Math.max(max, value);
            }
        }
        double range = Math.max(1e-9, max - min);
        int index = 0;
        for (double[] row : values) {
            for (double value : row) {
                pixels[index++] = (float) ((value - min) / range * 255.0);
            }
        }
        return new FloatProcessor(width, height, pixels);
    }

    private static double[][] crop(double[][] values, Box box) {
        double[][] out = new double[box.height()][box.width()];
        for (int y = 0; y < box.height(); y++) {
            System.arraycopy(values[box.top() + y], box.left(), out[y], 0, box.width());
        }
        return out;
    }

    private static boolean[][] makeDiskMask(int width, int height, double centerX, double centerY, double radius) {
        boolean[][] mask = new boolean[height][width];
        double radiusSquared = radius * radius;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double dx = x - centerX;
                double dy = y - centerY;
                mask[y][x] = dx * dx + dy * dy <= radiusSquared;
            }
        }
        return mask;
    }

    private static boolean[][] erodeDisk(boolean[][] mask, int radius) {
        if (radius <= 0) {
            return copy(mask);
        }
        boolean[][] out = new boolean[mask.length][mask[0].length];
        List<int[]> offsets = diskOffsets(radius);
        for (int y = 0; y < mask.length; y++) {
            for (int x = 0; x < mask[0].length; x++) {
                if (!mask[y][x]) {
                    continue;
                }
                boolean keep = true;
                for (int[] offset : offsets) {
                    int yy = y + offset[1];
                    int xx = x + offset[0];
                    if (yy < 0 || xx < 0 || yy >= mask.length || xx >= mask[0].length || !mask[yy][xx]) {
                        keep = false;
                        break;
                    }
                }
                out[y][x] = keep;
            }
        }
        return countTrue(out) > 0 ? out : copy(mask);
    }

    private static boolean[][] open(boolean[][] mask, int radius) {
        return dilate(erode(mask, radius), radius);
    }

    private static boolean[][] close(boolean[][] mask, int radius) {
        return erode(dilate(mask, radius), radius);
    }

    private static boolean[][] erode(boolean[][] mask, int radius) {
        boolean[][] out = new boolean[mask.length][mask[0].length];
        List<int[]> offsets = diskOffsets(radius);
        for (int y = 0; y < mask.length; y++) {
            for (int x = 0; x < mask[0].length; x++) {
                if (!mask[y][x]) {
                    continue;
                }
                boolean keep = true;
                for (int[] offset : offsets) {
                    int yy = y + offset[1];
                    int xx = x + offset[0];
                    if (yy < 0 || xx < 0 || yy >= mask.length || xx >= mask[0].length || !mask[yy][xx]) {
                        keep = false;
                        break;
                    }
                }
                out[y][x] = keep;
            }
        }
        return out;
    }

    private static boolean[][] dilate(boolean[][] mask, int radius) {
        boolean[][] out = new boolean[mask.length][mask[0].length];
        List<int[]> offsets = diskOffsets(radius);
        for (int y = 0; y < mask.length; y++) {
            for (int x = 0; x < mask[0].length; x++) {
                if (!mask[y][x]) {
                    continue;
                }
                for (int[] offset : offsets) {
                    int yy = y + offset[1];
                    int xx = x + offset[0];
                    if (yy >= 0 && xx >= 0 && yy < mask.length && xx < mask[0].length) {
                        out[yy][xx] = true;
                    }
                }
            }
        }
        return out;
    }

    private static List<int[]> diskOffsets(int radius) {
        List<int[]> offsets = new ArrayList<>();
        for (int y = -radius; y <= radius; y++) {
            for (int x = -radius; x <= radius; x++) {
                if (x * x + y * y <= radius * radius) {
                    offsets.add(new int[] {x, y});
                }
            }
        }
        return offsets;
    }

    private static boolean[][] fillHoles(boolean[][] mask) {
        int height = mask.length;
        int width = mask[0].length;
        boolean[][] visited = new boolean[height][width];
        ArrayDeque<int[]> queue = new ArrayDeque<>();
        for (int x = 0; x < width; x++) {
            enqueueIfHole(mask, visited, queue, 0, x);
            enqueueIfHole(mask, visited, queue, height - 1, x);
        }
        for (int y = 0; y < height; y++) {
            enqueueIfHole(mask, visited, queue, y, 0);
            enqueueIfHole(mask, visited, queue, y, width - 1);
        }
        int[][] neighbors = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
        while (!queue.isEmpty()) {
            int[] point = queue.removeFirst();
            for (int[] neighbor : neighbors) {
                int yy = point[0] + neighbor[0];
                int xx = point[1] + neighbor[1];
                if (yy < 0 || xx < 0 || yy >= height || xx >= width || visited[yy][xx] || mask[yy][xx]) {
                    continue;
                }
                visited[yy][xx] = true;
                queue.addLast(new int[] {yy, xx});
            }
        }
        boolean[][] filled = copy(mask);
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                if (!mask[y][x] && !visited[y][x]) {
                    filled[y][x] = true;
                }
            }
        }
        return filled;
    }

    private static void enqueueIfHole(boolean[][] mask, boolean[][] visited, ArrayDeque<int[]> queue, int y, int x) {
        if (!mask[y][x] && !visited[y][x]) {
            visited[y][x] = true;
            queue.addLast(new int[] {y, x});
        }
    }

    private static boolean[][] removeSmall(boolean[][] mask, int minArea) {
        boolean[][] out = new boolean[mask.length][mask[0].length];
        for (ConnectedComponent component : connectedComponents(mask)) {
            if (component.pixels().size() >= minArea) {
                for (int[] point : component.pixels()) {
                    out[point[0]][point[1]] = true;
                }
            }
        }
        return out;
    }

    private static List<ConnectedComponent> connectedComponents(boolean[][] mask) {
        int height = mask.length;
        int width = mask[0].length;
        boolean[][] visited = new boolean[height][width];
        List<ConnectedComponent> components = new ArrayList<>();
        int[][] neighbors = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                if (!mask[y][x] || visited[y][x]) {
                    continue;
                }
                ArrayDeque<int[]> queue = new ArrayDeque<>();
                List<int[]> pixels = new ArrayList<>();
                queue.add(new int[] {y, x});
                visited[y][x] = true;
                int minX = x;
                int minY = y;
                int maxX = x;
                int maxY = y;
                while (!queue.isEmpty()) {
                    int[] point = queue.removeFirst();
                    pixels.add(point);
                    minX = Math.min(minX, point[1]);
                    minY = Math.min(minY, point[0]);
                    maxX = Math.max(maxX, point[1]);
                    maxY = Math.max(maxY, point[0]);
                    for (int[] neighbor : neighbors) {
                        int yy = point[0] + neighbor[0];
                        int xx = point[1] + neighbor[1];
                        if (yy < 0 || xx < 0 || yy >= height || xx >= width || visited[yy][xx] || !mask[yy][xx]) {
                            continue;
                        }
                        visited[yy][xx] = true;
                        queue.addLast(new int[] {yy, xx});
                    }
                }
                components.add(new ConnectedComponent(pixels, new Rectangle(minX, minY, maxX - minX + 1, maxY - minY + 1)));
            }
        }
        return components;
    }

    private static double[][] distanceTransform(boolean[][] mask) {
        int height = mask.length;
        int width = mask[0].length;
        double[][] dist = new double[height][width];
        double inf = 1_000_000.0;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                dist[y][x] = mask[y][x] ? inf : 0.0;
            }
        }
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                if (!mask[y][x]) {
                    continue;
                }
                if (y > 0) {
                    dist[y][x] = Math.min(dist[y][x], dist[y - 1][x] + 1.0);
                    if (x > 0) {
                        dist[y][x] = Math.min(dist[y][x], dist[y - 1][x - 1] + Math.sqrt(2.0));
                    }
                    if (x + 1 < width) {
                        dist[y][x] = Math.min(dist[y][x], dist[y - 1][x + 1] + Math.sqrt(2.0));
                    }
                }
                if (x > 0) {
                    dist[y][x] = Math.min(dist[y][x], dist[y][x - 1] + 1.0);
                }
            }
        }
        for (int y = height - 1; y >= 0; y--) {
            for (int x = width - 1; x >= 0; x--) {
                if (!mask[y][x]) {
                    continue;
                }
                if (y + 1 < height) {
                    dist[y][x] = Math.min(dist[y][x], dist[y + 1][x] + 1.0);
                    if (x > 0) {
                        dist[y][x] = Math.min(dist[y][x], dist[y + 1][x - 1] + Math.sqrt(2.0));
                    }
                    if (x + 1 < width) {
                        dist[y][x] = Math.min(dist[y][x], dist[y + 1][x + 1] + Math.sqrt(2.0));
                    }
                }
                if (x + 1 < width) {
                    dist[y][x] = Math.min(dist[y][x], dist[y][x + 1] + 1.0);
                }
            }
        }
        return dist;
    }

    private static double[] collectValues(double[][] values, boolean[][] mask) {
        List<Double> out = new ArrayList<>();
        for (int y = 0; y < mask.length; y++) {
            for (int x = 0; x < mask[0].length; x++) {
                if (mask[y][x]) {
                    out.add(values[y][x]);
                }
            }
        }
        double[] array = out.stream().mapToDouble(Double::doubleValue).toArray();
        Arrays.sort(array);
        return array;
    }

    private static double quantile(double[] values, double q) {
        if (values.length == 0) {
            return 0.0;
        }
        double position = q * (values.length - 1);
        int lower = (int) Math.floor(position);
        int upper = (int) Math.ceil(position);
        if (lower == upper) {
            return values[lower];
        }
        double weight = position - lower;
        return values[lower] * (1.0 - weight) + values[upper] * weight;
    }

    private static double otsuThreshold(double[] values) {
        if (values.length == 0) {
            return 0.0;
        }
        double min = values[0];
        double max = values[values.length - 1];
        if (max <= min) {
            return min;
        }
        int[] histogram = new int[256];
        for (double value : values) {
            int index = clamp((int) Math.round((value - min) / (max - min) * 255.0), 0, 255);
            histogram[index]++;
        }
        AutoThresholder threshold = new AutoThresholder();
        int level = threshold.getThreshold(AutoThresholder.Method.Otsu, histogram);
        return min + (max - min) * (level / 255.0);
    }

    private static double[][] gaussian(double[][] image, double sigma) {
        int width = image[0].length;
        int height = image.length;
        float[] pixels = new float[width * height];
        int index = 0;
        for (double[] row : image) {
            for (double value : row) {
                pixels[index++] = (float) value;
            }
        }
        FloatProcessor processor = new FloatProcessor(width, height, pixels);
        new GaussianBlur().blurFloat(processor, sigma, sigma, 0.01);
        float[] blurred = (float[]) processor.getPixels();
        double[][] out = new double[height][width];
        index = 0;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                out[y][x] = blurred[index++];
            }
        }
        return out;
    }

    private static ByteProcessor toByteProcessor(boolean[][] mask) {
        int width = mask[0].length;
        int height = mask.length;
        byte[] pixels = new byte[width * height];
        int index = 0;
        for (boolean[] row : mask) {
            for (boolean value : row) {
                pixels[index++] = (byte) (value ? 255 : 0);
            }
        }
        return new ByteProcessor(width, height, pixels);
    }

    private static boolean[][] fromByteProcessor(ByteProcessor processor) {
        boolean[][] out = new boolean[processor.getHeight()][processor.getWidth()];
        byte[] pixels = (byte[]) processor.getPixels();
        int index = 0;
        for (int y = 0; y < processor.getHeight(); y++) {
            for (int x = 0; x < processor.getWidth(); x++) {
                out[y][x] = (pixels[index++] & 0xFF) != 0;
            }
        }
        return out;
    }

    private static boolean[][] copy(boolean[][] mask) {
        boolean[][] out = new boolean[mask.length][mask[0].length];
        for (int y = 0; y < mask.length; y++) {
            System.arraycopy(mask[y], 0, out[y], 0, mask[0].length);
        }
        return out;
    }

    private static int countTrue(boolean[][] mask) {
        int count = 0;
        for (boolean[] row : mask) {
            for (boolean value : row) {
                if (value) {
                    count++;
                }
            }
        }
        return count;
    }

    private static Metadata parseMetadata(Path imagePath, DatasetMode datasetMode) {
        String imageName = imagePath.getFileName().toString();
        String stem = stripExtension(imageName);
        int captureIndex = 1;
        String plateKey = stem;
        int matchIndex = stem.lastIndexOf('_');
        if (matchIndex >= 0) {
            String maybeIndex = stem.substring(matchIndex + 1);
            if (maybeIndex.chars().allMatch(Character::isDigit)) {
                captureIndex = Integer.parseInt(maybeIndex);
                plateKey = stem.substring(0, matchIndex);
            }
        }
        String[] tokens = plateKey.split("_");
        Double doseGy = tokens.length >= 1 ? parseDecimal(tokens[0]) : null;
        String doseUnit = tokens.length >= 2 ? tokens[1] : null;
        Double top = tokens.length >= 3 ? parseConcentrationToken(tokens[2]) : null;
        Double bottom = tokens.length >= 4 ? parseConcentrationToken(tokens[3]) : null;
        Map<Integer, Map<String, Object>> manifestRows = MANIFEST.rowsForImage(imageName);
        Map<String, Object> imageFields = new LinkedHashMap<>();
        if (!manifestRows.isEmpty()) {
            Map<String, Object> first = manifestRows.values().iterator().next();
            plateKey = stringField(first, "plate_key", plateKey);
            captureIndex = intField(first, "capture_index", captureIndex);
            doseGy = doubleField(first, "dose_gy", doseGy);
            doseUnit = stringField(first, "dose_unit", doseUnit);
            top = doubleField(first, "top_concentration", top);
            bottom = doubleField(first, "bottom_concentration", bottom);
            imageFields.putAll(collectCommonImageFields(manifestRows));
        }
        boolean reverseColumnAssignments = isZrDataset(imagePath, datasetMode) && !isZrNormalOrientationControl(doseGy, top, bottom);
        imageFields.put("image_name", imageName);
        imageFields.put("plate_key", plateKey);
        imageFields.put("capture_index", captureIndex);
        if (doseGy != null) {
            imageFields.put("dose_gy", formatNumberObject(doseGy));
        }
        if (doseUnit != null) {
            imageFields.put("dose_unit", doseUnit);
        }
        return new Metadata(imageName, plateKey, captureIndex, doseGy, doseUnit, top, bottom, reverseColumnAssignments, imageFields, manifestRows);
    }

    private static Double parseDecimal(String raw) {
        try {
            return Double.parseDouble(raw.replace(',', '.'));
        } catch (NumberFormatException ignored) {
            return null;
        }
    }

    private static Double parseConcentrationToken(String raw) {
        if (raw == null) {
            return null;
        }
        String normalized = raw.replace(',', '.');
        if ("0.75".equals(normalized)) {
            return 0.075;
        }
        return parseDecimal(raw);
    }

    private static boolean isZrDataset(Path imagePath, DatasetMode datasetMode) {
        if (datasetMode == DatasetMode.ZR) {
            return true;
        }
        if (datasetMode == DatasetMode.HF) {
            return false;
        }
        String fullPath = imagePath.toString().toLowerCase(Locale.ROOT);
        return fullPath.contains("clon zr egor") || fullPath.contains("/zr/") || fullPath.contains("\\zr\\");
    }

    private static boolean isZrNormalOrientationControl(Double doseGy, Double top, Double bottom) {
        return nearlyEquals(doseGy, 0.0) && nearlyEquals(top, 0.0) && nearlyEquals(bottom, 0.075);
    }

    private static boolean nearlyEquals(Double value, double expected) {
        return value != null && Math.abs(value - expected) < 1e-9;
    }

    private static Double concentrationForWell(Metadata metadata, int wellIndex) {
        boolean leftColumn = wellIndex <= 3;
        if (metadata.reverseColumnAssignments()) {
            return leftColumn ? metadata.bottomConcentration() : metadata.topConcentration();
        }
        return leftColumn ? metadata.topConcentration() : metadata.bottomConcentration();
    }

    private static WellAssignment resolveWellAssignment(Metadata metadata, int wellIndex) {
        String defaultRowName;
        Double defaultConcentration;
        int defaultReplicate;
        if (metadata.wellFields().isEmpty() && CUSTOM_LAYOUT.isEmpty()) {
            defaultRowName = wellIndex <= 3 ? "top" : "bottom";
            defaultConcentration = concentrationForWell(metadata, wellIndex);
            defaultReplicate = defaultRowName.equals("top") ? wellIndex : wellIndex - 3;
        } else {
            defaultRowName = "group_1";
            defaultConcentration = null;
            defaultReplicate = wellIndex;
        }

        Map<String, Object> manifestRow = metadata.wellFields().getOrDefault(wellIndex, Map.of());
        String rowName = stringField(manifestRow, "row_name", defaultRowName);
        Double concentration = doubleField(manifestRow, "concentration", defaultConcentration);
        int replicate = intField(manifestRow, "replicate", defaultReplicate);

        Map<String, Object> extraFields = new LinkedHashMap<>();
        extraFields.put("position_label", "W" + wellIndex);
        for (Map.Entry<String, Object> entry : manifestRow.entrySet()) {
            String key = entry.getKey();
            if (isImageLevelField(key) || isStandardWellField(key)) {
                continue;
            }
            extraFields.put(key, entry.getValue());
        }
        return new WellAssignment(rowName, concentration, replicate, extraFields);
    }

    private static boolean isStandardWellField(String key) {
        return "image_name".equals(key)
            || "well_index".equals(key)
            || "row_name".equals(key)
            || "replicate".equals(key)
            || "concentration".equals(key);
    }

    private static boolean isImageLevelField(String key) {
        return "plate_key".equals(key)
            || "capture_index".equals(key)
            || "dose_gy".equals(key)
            || "dose_unit".equals(key)
            || "top_concentration".equals(key)
            || "bottom_concentration".equals(key);
    }

    private static Object formatNumberObject(Double value) {
        if (value == null) {
            return null;
        }
        return round(value, 6);
    }

    private static Map<String, Object> collectCommonImageFields(Map<Integer, Map<String, Object>> manifestRows) {
        if (manifestRows.isEmpty()) {
            return Map.of();
        }
        Map<String, Object> first = manifestRows.values().iterator().next();
        Map<String, Object> common = new LinkedHashMap<>();
        for (Map.Entry<String, Object> entry : first.entrySet()) {
            String key = entry.getKey();
            if ("well_index".equals(key) || "image_name".equals(key)) {
                continue;
            }
            boolean same = true;
            for (Map<String, Object> row : manifestRows.values()) {
                Object other = row.get(key);
                if (other == null || !other.toString().equals(entry.getValue().toString())) {
                    same = false;
                    break;
                }
            }
            if (same) {
                common.put(key, entry.getValue());
            }
        }
        return common;
    }

    private static String stringField(Map<String, Object> row, String key, String defaultValue) {
        Object value = row.get(key);
        if (value == null || value.toString().isBlank()) {
            return defaultValue;
        }
        return value.toString().trim();
    }

    private static int intField(Map<String, Object> row, String key, int defaultValue) {
        Object value = row.get(key);
        if (value == null || value.toString().isBlank()) {
            return defaultValue;
        }
        return Integer.parseInt(value.toString().trim());
    }

    private static Double doubleField(Map<String, Object> row, String key, Double defaultValue) {
        Object value = row.get(key);
        if (value == null || value.toString().isBlank()) {
            return defaultValue;
        }
        return parseDecimal(value.toString());
    }

    private static String stripExtension(String name) {
        int index = name.lastIndexOf('.');
        return index >= 0 ? name.substring(0, index) : name;
    }

    private static Map<String, Object> baseMetadataMap(Metadata metadata) {
        Map<String, Object> map = new LinkedHashMap<>();
        map.put("image_name", metadata.imageName());
        map.put("plate_key", metadata.plateKey());
        map.put("capture_index", metadata.captureIndex());
        map.put("dose_gy", metadata.doseGy());
        map.put("dose_unit", metadata.doseUnit());
        map.put("top_concentration", metadata.topConcentration());
        map.put("bottom_concentration", metadata.bottomConcentration());
        map.put("reverse_column_assignments", metadata.reverseColumnAssignments());
        for (Map.Entry<String, Object> entry : metadata.imageFields().entrySet()) {
            map.put(entry.getKey(), entry.getValue());
        }
        return map;
    }

    private static Map<String, Object> metadataToMap(Metadata metadata) {
        return baseMetadataMap(metadata);
    }

    private static List<Map<String, Object>> wellSummariesToMaps(List<WellSummary> summaries) {
        List<Map<String, Object>> rows = new ArrayList<>();
        for (WellSummary summary : summaries) {
            Map<String, Object> row = new LinkedHashMap<>();
            row.putAll(summary.extraFields());
            row.put("well_index", summary.wellIndex());
            row.put("row_name", summary.rowName());
            row.put("replicate", summary.replicate());
            row.put("concentration", summary.concentration());
            row.put("raw_component_count", summary.rawComponentCount());
            row.put("counted_colonies", summary.countedColonies());
            row.put("single_area_px", round(summary.singleAreaPx(), 3));
            row.put("threshold", round(summary.threshold(), 6));
            row.put("well_area_px", summary.wellAreaPx());
            row.put("count_area_px", summary.countAreaPx());
            rows.add(row);
        }
        return rows;
    }

    private static void writeCsv(Path path, List<Map<String, Object>> rows) throws IOException {
        if (rows.isEmpty()) {
            return;
        }
        List<String> columns = new ArrayList<>(rows.getFirst().keySet());
        columns.sort(Comparator.naturalOrder());
        try (BufferedWriter writer = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            writer.write(String.join(",", columns));
            writer.write("\n");
            for (Map<String, Object> row : rows) {
                for (int i = 0; i < columns.size(); i++) {
                    if (i > 0) {
                        writer.write(",");
                    }
                    writer.write(csvValue(row.get(columns.get(i))));
                }
                writer.write("\n");
            }
        }
    }

    private static String csvValue(Object value) {
        if (value == null) {
            return "";
        }
        String raw = value.toString();
        if (raw.contains(",") || raw.contains("\"")) {
            return "\"" + raw.replace("\"", "\"\"") + "\"";
        }
        return raw;
    }

    private static void writeGroupSummary(List<Map<String, Object>> wellRows, Path path) throws IOException {
        List<String> groupFields = resolveSummaryGroupFields(wellRows);
        Map<String, List<Integer>> groups = new HashMap<>();
        Map<String, Map<String, Object>> meta = new HashMap<>();
        for (Map<String, Object> row : wellRows) {
            String key = buildKey(row, groupFields);
            groups.computeIfAbsent(key, unused -> new ArrayList<>()).add(((Number) row.get("counted_colonies")).intValue());
            meta.putIfAbsent(key, row);
        }
        List<Map<String, Object>> out = new ArrayList<>();
        for (Map.Entry<String, List<Integer>> entry : groups.entrySet()) {
            List<Integer> counts = entry.getValue();
            Map<String, Object> source = meta.get(entry.getKey());
            double mean = counts.stream().mapToInt(Integer::intValue).average().orElse(0.0);
            double variance = 0.0;
            for (int count : counts) {
                variance += Math.pow(count - mean, 2.0);
            }
            variance /= Math.max(1, counts.size());
            Map<String, Object> row = new LinkedHashMap<>();
            for (String field : groupFields) {
                row.put(field, source.get(field));
            }
            row.put("replicates", counts.size());
            row.put("mean_count", round(mean, 3));
            row.put("std_count", round(Math.sqrt(variance), 3));
            row.put("counts", counts.toString().replace(" ", ""));
            out.add(row);
        }
        writeCsv(path, out);
    }

    private static void writeDuplicateWellSummary(List<Map<String, Object>> wellRows, Path path) throws IOException {
        Map<String, List<Map<String, Object>>> groups = new HashMap<>();
        for (Map<String, Object> row : wellRows) {
            String key = row.get("plate_key") + "|" + row.get("well_index");
            groups.computeIfAbsent(key, unused -> new ArrayList<>()).add(row);
        }
        List<Map<String, Object>> out = new ArrayList<>();
        for (List<Map<String, Object>> rows : groups.values()) {
            rows.sort(Comparator.comparingInt(row -> ((Number) row.get("capture_index")).intValue()));
            double[] counts = rows.stream().mapToDouble(row -> ((Number) row.get("counted_colonies")).doubleValue()).toArray();
            double mean = Arrays.stream(counts).average().orElse(0.0);
            double variance = 0.0;
            for (double count : counts) {
                variance += Math.pow(count - mean, 2.0);
            }
            variance /= Math.max(1, counts.length);
            Map<String, Object> first = rows.getFirst();
            Map<String, Object> row = new HashMap<>();
            row.put("plate_key", first.get("plate_key"));
            row.put("dose_gy", first.get("dose_gy"));
            row.put("dose_unit", first.get("dose_unit"));
            row.put("top_concentration", first.get("top_concentration"));
            row.put("bottom_concentration", first.get("bottom_concentration"));
            row.put("well_index", first.get("well_index"));
            row.put("row_name", first.get("row_name"));
            row.put("replicate", first.get("replicate"));
            row.put("concentration", first.get("concentration"));
            row.put("image_count", rows.size());
            row.put("mean_count", round(mean, 3));
            row.put("std_count", round(Math.sqrt(variance), 3));
            row.put("min_count", Arrays.stream(counts).min().orElse(0.0));
            row.put("max_count", Arrays.stream(counts).max().orElse(0.0));
            row.put("counts", rows.stream().map(value -> value.get("counted_colonies").toString()).toList().toString().replace(" ", ""));
            row.put("capture_indices", rows.stream().map(value -> value.get("capture_index").toString()).toList().toString().replace(" ", ""));
            row.put("images", rows.stream().map(value -> value.get("image_name").toString()).toList().toString().replace(" ", ""));
            out.add(row);
        }
        writeCsv(path, out);
    }

    private static void writeDuplicateGroupSummary(List<Map<String, Object>> wellRows, Path path) throws IOException {
        List<String> groupFields = resolveSummaryGroupFields(wellRows);
        Map<String, List<Map<String, Object>>> groups = new HashMap<>();
        for (Map<String, Object> row : wellRows) {
            String key = row.get("plate_key") + "|" + buildKey(row, groupFields);
            groups.computeIfAbsent(key, unused -> new ArrayList<>()).add(row);
        }
        List<Map<String, Object>> out = new ArrayList<>();
        for (List<Map<String, Object>> rows : groups.values()) {
            double[] counts = rows.stream().mapToDouble(row -> ((Number) row.get("counted_colonies")).doubleValue()).toArray();
            double mean = Arrays.stream(counts).average().orElse(0.0);
            double variance = 0.0;
            for (double count : counts) {
                variance += Math.pow(count - mean, 2.0);
            }
            variance /= Math.max(1, counts.length);
            Map<String, Object> first = rows.getFirst();
            Set<String> images = new HashSet<>();
            for (Map<String, Object> row : rows) {
                images.add(row.get("image_name").toString());
            }
            Map<String, Object> outRow = new LinkedHashMap<>();
            outRow.put("plate_key", first.get("plate_key"));
            outRow.put("dose_gy", first.get("dose_gy"));
            outRow.put("dose_unit", first.get("dose_unit"));
            for (String field : groupFields) {
                outRow.put(field, first.get(field));
            }
            outRow.put("image_count", images.size());
            outRow.put("replicate_wells", rows.size());
            outRow.put("mean_count", round(mean, 3));
            outRow.put("std_count", round(Math.sqrt(variance), 3));
            outRow.put("counts", rows.stream().map(value -> value.get("counted_colonies").toString()).toList().toString().replace(" ", ""));
            outRow.put("images", images.toString().replace(" ", ""));
            out.add(outRow);
        }
        writeCsv(path, out);
    }

    private static List<String> resolveSummaryGroupFields(List<Map<String, Object>> rows) {
        if (rows.isEmpty()) {
            return List.of();
        }
        if (!SUMMARY_GROUP_FIELDS.isEmpty()) {
            List<String> resolved = new ArrayList<>();
            for (String field : SUMMARY_GROUP_FIELDS) {
                if (rows.getFirst().containsKey(field)) {
                    resolved.add(field);
                }
            }
            if (!resolved.isEmpty()) {
                return resolved;
            }
        }
        Set<String> excluded = Set.of(
            "image_name",
            "plate_key",
            "capture_index",
            "top_concentration",
            "bottom_concentration",
            "reverse_column_assignments",
            "well_index",
            "replicate",
            "row_name",
            "position_label",
            "raw_component_count",
            "counted_colonies",
            "single_area_px",
            "threshold",
            "well_area_px",
            "count_area_px",
            "analysis_source",
            "analysis_offset_x_px",
            "analysis_offset_y_px",
            "analysis_scale",
            "well_center_x_px",
            "well_center_y_px",
            "global_well_center_x_px",
            "global_well_center_y_px",
            "well_radius_px",
            "rough_center_x_px",
            "rough_center_y_px",
            "global_rough_center_x_px",
            "global_rough_center_y_px",
            "rough_radius_px"
        );
        List<String> inferred = new ArrayList<>();
        for (String key : rows.getFirst().keySet()) {
            if (!excluded.contains(key)) {
                inferred.add(key);
            }
        }
        if (inferred.isEmpty()) {
            return List.of("well_index");
        }
        return inferred;
    }

    private static String buildKey(Map<String, Object> row, List<String> fields) {
        StringBuilder builder = new StringBuilder();
        for (String field : fields) {
            builder.append(field).append('=').append(row.getOrDefault(field, "")).append('|');
        }
        return builder.toString();
    }

    private static void writeJson(Path path, Map<String, Object> data) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            writer.write(toJson(data, 0));
            writer.write("\n");
        }
    }

    private static String toJson(Object value, int indent) {
        String pad = " ".repeat(indent);
        if (value == null) {
            return "null";
        }
        if (value instanceof String string) {
            return "\"" + string.replace("\\", "\\\\").replace("\"", "\\\"") + "\"";
        }
        if (value instanceof Number || value instanceof Boolean) {
            return value.toString();
        }
        if (value instanceof Map<?, ?> map) {
            StringBuilder builder = new StringBuilder();
            builder.append("{\n");
            int index = 0;
            for (Map.Entry<?, ?> entry : map.entrySet()) {
                builder.append(" ".repeat(indent + 2))
                    .append(toJson(entry.getKey().toString(), 0))
                    .append(": ")
                    .append(toJson(entry.getValue(), indent + 2));
                if (++index < map.size()) {
                    builder.append(",");
                }
                builder.append("\n");
            }
            builder.append(pad).append("}");
            return builder.toString();
        }
        if (value instanceof List<?> list) {
            StringBuilder builder = new StringBuilder();
            builder.append("[");
            for (int i = 0; i < list.size(); i++) {
                if (i > 0) {
                    builder.append(", ");
                }
                builder.append(toJson(list.get(i), indent + 2));
            }
            builder.append("]");
            return builder.toString();
        }
        return toJson(value.toString(), indent);
    }

    private static double round(double value, int places) {
        double factor = Math.pow(10.0, places);
        return Math.round(value * factor) / factor;
    }

    private static int clamp(int value, int min, int max) {
        return Math.max(min, Math.min(max, value));
    }
}
