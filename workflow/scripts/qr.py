# Copyright 2026 Xin Huang and Andrea Koca
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import pandas as pd
from quantile_forest import RandomForestQuantileRegressor
from sklearn.linear_model import QuantileRegressor
from sklearn.ensemble import GradientBoostingRegressor


def get_xy(df: pd.DataFrame, features: list[str]) -> tuple[pd.DataFrame, pd.Series]:
    """
    Extract model features and target.
    """
    x = df[features].copy()
    y = df["S*_score"].copy()
    return x, y


def get_model(model_type: str, quantile: float):
    """
    Return the selected regression model.
    """
    if model_type == "quantile":

        return QuantileRegressor(
            quantile=quantile,
            alpha=0,
            solver="highs",
        )

    if model_type == "gradient":

        return GradientBoostingRegressor(
            loss="quantile",
            alpha=quantile,
            n_estimators=200,
            max_depth=3,
        )

    if model_type == "qrf":

        return RandomForestQuantileRegressor(
            n_estimators=200,
            n_jobs=-1,
        )

    raise ValueError(
        f"Unsupported model_type: {model_type}. "
        "Choose from: 'quantile', 'gradient', or 'qrf'."
    )


quantile = float(snakemake.wildcards.quantile)
model_type = snakemake.wildcards.qr_model
feature_set = snakemake.wildcards.feature_set


null_df = pd.read_csv(snakemake.input.training_data, sep="\t")
score_df = pd.read_csv(snakemake.input.test_data, sep="\t")

feature_sets = {
    "sstar_snp": ["S*_SNP_number"],
    "region_snp": ["region_ind_SNP_number"],
    "both": ["S*_SNP_number", "region_ind_SNP_number"],
}


if feature_set not in feature_sets:
    raise ValueError(
        f"Unsupported feature_set: {feature_set}. "
        f"Choose from: {list(feature_sets.keys())}"
    )

features = feature_sets[feature_set]

train = null_df.dropna(subset=features + ["S*_score"]).copy()
pred_df = score_df.dropna(subset=features + ["S*_score"]).copy()

x_train, y_train = get_xy(train, features)
x_pred, _ = get_xy(pred_df, features)

model = get_model(model_type, quantile)
model.fit(x_train, y_train)

pred_df = pred_df.copy()

if model_type == "qrf":
    pred_df["expected_S*_score"] = model.predict(
        x_pred,
        quantiles=quantile,
    )
else:
    pred_df["expected_S*_score"] = model.predict(x_pred)


pred_df.to_csv(snakemake.output.preds, sep="\t", index=False)
