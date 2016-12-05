<?php

/**
 * Vincenty formula to calculates distance between two points.
 */
function gpx_parser_calculate_vincenty($lat1, $lon1, $lat2, $lon2) {
  $lat1 = deg2rad($lat1);
  $lat2 = deg2rad($lat2);
  $lon1 = deg2rad($lon1);
  $lon2 = deg2rad($lon2);
  $a = 6378137;
  $b = 6356752.314245;
  $f = 1 / 298.257223563;
  // WGS-84 ellipsoid.
  $londiff = $lon2 - $lon1;
  $u1 = atan((1 - $f) * tan($lat1));
  $u2 = atan((1 - $f) * tan($lat2));
  $sin_u1 = sin($u1);
  $cos_u1 = cos($u1);
  $sin_u2 = sin($u2);
  $cos_u2 = cos($u2);
  $lambda = $londiff;
  $lambda_p = 2 * M_PI;
  $iter_limit = 100;
  while (abs($lambda - $lambda_p) > 1e-12 && --$iter_limit > 0) {
    $sin_lambda = sin($lambda);
    $cos_lambda = cos($lambda);
    $sin_sigma = sqrt(($cos_u2 * $sin_lambda) * ($cos_u2 * $sin_lambda) + ($cos_u1 * $sin_u2 - $sin_u1 * $cos_u2 * $cos_lambda) * ($cos_u1 * $sin_u2 - $sin_u1 * $cos_u2 * $cos_lambda));
    if ($sin_sigma == 0) {
      return 0;
    }
    // Co-incident points.
    $cos_sigma = $sin_u1 * $sin_u2 + $cos_u1 * $cos_u2 * $cos_lambda;
    $sigma = atan2($sin_sigma, $cos_sigma);
    // Was atan2.
    $alpha = asin($cos_u1 * $cos_u2 * $sin_lambda / $sin_sigma);
    $cos_sq_alpha = cos($alpha) * cos($alpha);
    $cos2sigma_m = $cos_sigma - 2 * $sin_u1 * $sin_u2 / $cos_sq_alpha;
    $cc = $f / 16 * $cos_sq_alpha * (4 + $f * (4 - 3 * $cos_sq_alpha));
    $lambda_p = $lambda;
    $lambda = $londiff + (1 - $cc) * $f * sin($alpha) * ($sigma + $cc * $sin_sigma * ($cos2sigma_m + $cc * $cos_sigma * (-1 + 2 * $cos2sigma_m * $cos2sigma_m)));
  }
  if ($iter_limit == 0) {
    return FALSE;
  }
  // Formula failed to converge.
  $u_sq = $cos_sq_alpha * ($a * $a - $b * $b) / ($b * $b);
  $aa   = 1 + $u_sq / 16384 * (4096 + $u_sq * (-768 + $u_sq * (320 - 175 * $u_sq)));
  $bb   = $u_sq / 1024 * (256 + $u_sq * (-128 + $u_sq * (74 - 47 * $u_sq)));
  $delta_sigma = $bb * $sin_sigma * ($cos2sigma_m + $bb / 4 * ($cos_sigma * (-1 + 2 * $cos2sigma_m * $cos2sigma_m) -
        $bb / 6 * $cos2sigma_m * (-3 + 4 * $sin_sigma * $sin_sigma) * (-3 + 4 * $cos2sigma_m * $cos2sigma_m)));
  $s = $b * $aa * ($sigma - $delta_sigma);
  $s = round($s, 3);
  // Round to 1mm precision.
  return $s;
}
