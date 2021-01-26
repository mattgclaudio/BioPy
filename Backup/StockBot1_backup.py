#!/usr/bin/python3
from yahoo_fin.stock_info import *
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import numpy as np
import tensorflow as tf
from tensorflow import keras
from keras.models import Sequential
from keras.layers import Dense, LSTM
from datetime import date, datetime
import csv
import math
import sys

# GPU config
physical_devices = tf.config.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(physical_devices[0], enable=True)

scl = MinMaxScaler()

plt.style.use('bmh')
pd.options.display.width = 0


def predict_price(ticker):
    try:
        all_data = get_data(ticker, start_date="2010-1-1", end_date=date.today()).filter(['close'])
        # Numpy Array, no dates!
        close_set = all_data.values

        if np.isnan(close_set).any():
            print("we got naaans in here:\t" + ticker)
            # TODO only strip naan portions if there are gaps.

        else:
            scaled_set = scl.fit_transform(close_set)
            trainlength = math.ceil((len(close_set) * .8))
            trainset = scaled_set[0:trainlength, :]

            # X train: Closing costs 60 days prior to day D0 trainset[j - 60:j, 0] |||| Y Train: closing cost for day j
            x_train = []
            y_train = []

            for j in range(60, len(trainset)):
                x_train.append(trainset[j - 60:j, 0])
                y_train.append(trainset[j, 0])

            # list -> np.array
            x_train, y_train = np.array(x_train), np.array(y_train)

            # LSTM model requires 3D tensor with  shape [batch, timestep, features]
            x_train = np.reshape(x_train, (x_train.shape[0], x_train.shape[1], 1))

            # Building a new model.
            # close_model = Sequential()
            # close_model.add(LSTM(units=70, return_sequences=True))
            # close_model.add(LSTM(units=70, return_sequences=False))
            # close_model.add(Dense(units=35))
            # close_model.add(Dense(units=1))

            # Load trained model
            close_model = keras.models.load_model('model_docker')

            close_model.compile(optimizer=tf.keras.optimizers.SGD(learning_rate=0.01),
                                loss=tf.keras.losses.MeanSquaredError(),
                                metrics=['mse']
                                )
            close_model.fit(x_train, y_train, batch_size=1, epochs=1)

            close_model.save("model_docker")

            testset = scaled_set[trainlength - 60:, :]
            x_test = []
            y_test = close_set[trainlength:, :]  # valid closing cost test data

            for q in range(60, len(testset)):
                x_test.append(testset[q - 60:q, 0])  # 60 days prior --> present day exclusive

            x_test = np.array(x_test)
            x_test = np.reshape(x_test, (x_test.shape[0], x_test.shape[1], 1))

            guess = close_model.predict(x_test)
            inv_guess = scl.inverse_transform(guess)

            t = all_data[:trainlength]
            v = all_data[trainlength:]

            v['Predictions'] = scl.inverse_transform(guess)

            plt.figure(figsize=(10, 5))
            plt.title(ticker + ' Predictions')
            plt.xlabel('Date')
            plt.ylabel('Stock Price', fontsize=18)
            plt.plot(t['close'])
            plt.plot(v[['close', 'Predictions']])
            plt.legend(['Train', 'Valid', 'Predictions'])

            now = datetime.now()  # tag for graph
            figure_tag = now.strftime("%H:%M:%S") + "_" + ticker + "_Graph"
            plt.savefig(figure_tag)
            rmse_error = np.sqrt(np.mean(((inv_guess - y_test) ** 2)))
            print(rmse_error)

    except AssertionError:  # if the function above returns no data
        print("\n\nNo stock data for :\t" + ticker + "\n\n")

# symbols.csv len 471
symbol_file = pd.read_csv('symbols.csv')
stock_tickers = symbol_file.values[:,1]


for i in range(135, 139):
    print(predict_price(stock_tickers[i][1]))
