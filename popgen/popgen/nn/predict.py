import tensorflow as tf
from tensorflow import keras

def predict_keras(features, model_path):
    nn_model = keras.models.load_model(model_path)
    preds = nn_model.predict(features)
    return preds

