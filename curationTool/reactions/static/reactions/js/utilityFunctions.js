// Generates a list of distinct colors in HSL format, avoiding yellow hues, based on the specified count.
function getDistinctColors(count) {
    var colors = [];
    var hueStep = 360 / count;
    var lightness = 40; // Adjust for darker colors, but not too dark

    for (var i = 0; i < 360; i += hueStep) {
        // Skip yellow and near-yellow hues
        if (i >= 50 && i <= 70) {
            continue;
        }

        colors.push('hsl(' + i + ', 100%, ' + lightness + '%)');
    }
    return colors;
}
