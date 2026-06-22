/* USER CODE BEGIN Header */
/**
 ******************************************************************************
 * @file           : main.c
 * @brief          : Main program body
 ******************************************************************************
 * @attention
 *
 * Copyright (c) 2026 STMicroelectronics.
 * All rights reserved.
 *
 * This software is licensed under terms that can be found in the LICENSE file
 * in the root directory of this software component.
 * If no LICENSE file comes with this software, it is provided AS-IS.
 *
 ******************************************************************************
 */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "adc.h"
#include "dac.h"
#include "dma.h"
#include "tim.h"
#include "gpio.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ola.h"
#include "arm_math.h"
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */
#define FIR_NUM_TAPS 80
#define POLYPHASE_TAPS 20
#define UPSAMPLE_FACTOR 4
#define DOWNSAMPLE_FACTOR 4

typedef struct _Biquad
{
  float b0, b1, b2;
  float a1, a2;
  float x1, x2;
  float y1, y2;
} Biquad;
/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
#define BUFFER_SIZE 256
#define UPSAMPLE_BUFFER_SIZE (UPSAMPLE_FACTOR * BUFFER_SIZE)
// FIR_NUM_TAPS defined in the private typedef section; it is used by FIR struct

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */
Biquad biquad = {
    .b0 = 0.00391607f,
    .b1 = 0.00783215f,
    .b2 = 0.00391607f,
    .a1 = -1.81531792f,
    .a2 = 0.83098222f,
    .x1 = 0.0f,
    .x2 = 0.0f,
    .y1 = 0.0f,
    .y2 = 0.0f};

static float decim_state[FIR_NUM_TAPS + UPSAMPLE_BUFFER_SIZE - 1];
static float interp_state[FIR_NUM_TAPS / UPSAMPLE_FACTOR + BUFFER_SIZE - 1];

const float h_decimation_coeffs[FIR_NUM_TAPS] = {
    -0.00046170067f,
    -0.0006826209f,
    -0.00052289446f,
    9.0989922e-19f,
    0.00066670505f,
    0.0010935766f,
    0.00090502407f,
    -1.4736314e-18f,
    -0.001250687f,
    -0.0020762829f,
    -0.0017179494f,
    2.4273441e-18f,
    0.0023232347f,
    0.0037952333f,
    0.0030862971f,
    -3.6728746e-18f,
    -0.0040320467f,
    -0.0064817758f,
    -0.0051934653f,
    5.0820239e-18f,
    0.0066164049f,
    0.010529813f,
    0.0083678022f,
    -6.5097524e-18f,
    -0.010550481f,
    -0.016761257f,
    -0.013331607f,
    7.8091081e-18f,
    0.017002164f,
    0.02733178f,
    0.022115644f,
    -8.846352e-18f,
    -0.029864374f,
    -0.050256344f,
    -0.043358808f,
    9.5147237e-18f,
    0.0740235f,
    0.15820688f,
    0.22474334f,
    0.25f,
    0.22474334f,
    0.15820688f,
    0.0740235f,
    9.5147237e-18f,
    -0.043358808f,
    -0.050256344f,
    -0.029864374f,
    -8.846352e-18f,
    0.022115644f,
    0.02733178f,
    0.017002164f,
    7.8091081e-18f,
    -0.013331607f,
    -0.016761257f,
    -0.010550481f,
    -6.5097524e-18f,
    0.0083678022f,
    0.010529813f,
    0.0066164049f,
    5.0820239e-18f,
    -0.0051934653f,
    -0.0064817758f,
    -0.0040320467f,
    -3.6728746e-18f,
    0.0030862971f,
    0.0037952333f,
    0.0023232347f,
    2.4273441e-18f,
    -0.0017179494f,
    -0.0020762829f,
    -0.001250687f,
    -1.4736314e-18f,
    0.00090502407f,
    0.0010935766f,
    0.00066670505f,
    9.0989922e-19f,
    -0.00052289446f,
    -0.0006826209f,
    -0.00046170067f,
    0.0f};
const float h_interpolation_coeffs[FIR_NUM_TAPS] = {
    -0.0018468027f, -0.0027304836f, -0.0020915779f, 3.6395969e-18f,
    0.0026668202f, 0.0043743064f, 0.0036200963f, -5.8945254e-18f,
    -0.0050027479f, -0.0083051317f, -0.0068717977f, 9.7093765e-18f,
    0.0092929389f, 0.015180933f, 0.012345189f, -1.4691498e-17f,
    -0.016128187f, -0.025927103f, -0.020773861f, 2.0328096e-17f,
    0.026465619f, 0.04211925f, 0.033471209f, -2.603901e-17f,
    -0.042201922f, -0.067045027f, -0.053326428f, 3.1236432e-17f,
    0.068008656f, 0.10932712f, 0.088462575f, -3.5385408e-17f,
    -0.1194575f, -0.20102538f, -0.17343523f, 3.8058895e-17f,
    0.296094f, 0.63282751f, 0.89897337f, 1.0f,
    0.89897337f, 0.63282751f, 0.296094f, 3.8058895e-17f,
    -0.17343523f, -0.20102538f, -0.1194575f, -3.5385408e-17f,
    0.088462575f, 0.10932712f, 0.068008656f, 3.1236432e-17f,
    -0.053326428f, -0.067045027f, -0.042201922f, -2.603901e-17f,
    0.033471209f, 0.04211925f, 0.026465619f, 2.0328096e-17f,
    -0.020773861f, -0.025927103f, -0.016128187f, -1.4691498e-17f,
    0.012345189f, 0.015180933f, 0.0092929389f, 9.7093765e-18f,
    -0.0068717977f, -0.0083051317f, -0.0050027479f, -5.8945254e-18f,
    0.0036200963f, 0.0043743064f, 0.0026668202f, 3.6395969e-18f,
    -0.0020915779f, -0.0027304836f, -0.0018468027f, 0.0f};
OLA_State s;
float processInputBuffer[BUFFER_SIZE] = {0.0f};
float processOutputBuffer[BUFFER_SIZE] = {0.0f};

static arm_fir_decimate_instance_f32 decim;
static arm_fir_interpolate_instance_f32 interp;

float oversampledInputBuffer[UPSAMPLE_BUFFER_SIZE] = {0.0f};
float oversampledOutputBuffer[UPSAMPLE_BUFFER_SIZE] = {0.0f};
volatile uint16_t adc_buffer[UPSAMPLE_BUFFER_SIZE * 2] = {0};
volatile uint16_t out_buffer[UPSAMPLE_BUFFER_SIZE * 2] = {0};
volatile uint16_t pot_params[3] = {0};

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
/* USER CODE BEGIN PFP */
void process_block();
/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

/* USER CODE END 0 */

/**
 * @brief  The application entry point.
 * @retval int
 */
int main(void)
{

  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DMA_Init();
  MX_DAC_Init();
  MX_ADC1_Init();
  MX_TIM2_Init();
  MX_ADC2_Init();
  MX_TIM3_Init();
  MX_TIM6_Init();
  /* USER CODE BEGIN 2 */
  OLA_Init(&s);
  if (arm_fir_decimate_init_f32(&decim, FIR_NUM_TAPS, DOWNSAMPLE_FACTOR,
                                h_decimation_coeffs, decim_state, UPSAMPLE_BUFFER_SIZE) != ARM_MATH_SUCCESS)
    Error_Handler();
  if (arm_fir_interpolate_init_f32(&interp, UPSAMPLE_FACTOR, FIR_NUM_TAPS,
                                   h_interpolation_coeffs, interp_state, BUFFER_SIZE) != ARM_MATH_SUCCESS)
    Error_Handler();

  HAL_ADC_Start_DMA(&hadc2, (uint32_t *)pot_params, 3);
  HAL_ADC_Start_DMA(&hadc1, (uint32_t *)adc_buffer, UPSAMPLE_BUFFER_SIZE * 2);
  HAL_DAC_Start_DMA(&hdac, DAC_CHANNEL_1,
                    (uint32_t *)out_buffer, UPSAMPLE_BUFFER_SIZE * 2, DAC_ALIGN_12B_R);

  HAL_TIM_Base_Start(&htim2); // ADC timer
  HAL_TIM_Base_Start(&htim3); // Potentiometer polling timer

  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
 * @brief System Clock Configuration
 * @retval None
 */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
   */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE1);

  /** Initializes the RCC Oscillators according to the specified parameters
   * in the RCC_OscInitTypeDef structure.
   */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
  RCC_OscInitStruct.HSEState = RCC_HSE_ON;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSE;
  RCC_OscInitStruct.PLL.PLLM = 4;
  RCC_OscInitStruct.PLL.PLLN = 180;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV2;
  RCC_OscInitStruct.PLL.PLLQ = 2;
  RCC_OscInitStruct.PLL.PLLR = 2;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Activate the Over-Drive mode
   */
  if (HAL_PWREx_EnableOverDrive() != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
   */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK | RCC_CLOCKTYPE_SYSCLK | RCC_CLOCKTYPE_PCLK1 | RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV4;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV2;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_5) != HAL_OK)
  {
    Error_Handler();
  }
}

/* USER CODE BEGIN 4 */

void HAL_ADC_ConvHalfCpltCallback(ADC_HandleTypeDef *hadc)
{
  if (hadc == &hadc1)
  {
    process_block(&adc_buffer[0], &out_buffer[0]);
  }
}
void HAL_ADC_ConvCpltCallback(ADC_HandleTypeDef *hadc)
{
  if (hadc == &hadc1)
  {
    process_block(&adc_buffer[UPSAMPLE_BUFFER_SIZE], &out_buffer[UPSAMPLE_BUFFER_SIZE]);
  }
  else if (hadc == &hadc2)
  {
    HAL_ADC_Start_DMA(&hadc2, (uint32_t *)pot_params, 3);
  }
}

void process_block(uint16_t *input_buffer, uint16_t *output_buffer)
{

  HAL_GPIO_WritePin(GPIOA, GPIO_PIN_0, GPIO_PIN_SET);

  for (int i = 0; i < UPSAMPLE_BUFFER_SIZE; ++i)
  {
    oversampledInputBuffer[i] = (input_buffer[i] - 2048) / 2048.0f;
  }

  HAL_GPIO_WritePin(GPIOA, GPIO_PIN_2, GPIO_PIN_SET);
  arm_fir_decimate_f32(&decim, oversampledInputBuffer, processInputBuffer, UPSAMPLE_BUFFER_SIZE);
  HAL_GPIO_WritePin(GPIOA, GPIO_PIN_2, GPIO_PIN_RESET);

 HAL_GPIO_WritePin(GPIOA, GPIO_PIN_3, GPIO_PIN_SET);
  OLA_Process(&s, processInputBuffer, processOutputBuffer, BUFFER_SIZE);
  HAL_GPIO_WritePin(GPIOA, GPIO_PIN_3, GPIO_PIN_RESET);
  /* for(int i = 0; i < BUFFER_SIZE; ++i)
  {
    processOutputBuffer[i] =
        biquad.b0 * processInputBuffer[i]
      + biquad.b1 * biquad.x1
      + biquad.b2 * biquad.x2
      - biquad.a1 * biquad.y1
      - biquad.a2 * biquad.y2;

    biquad.x2 = biquad.x1;
    biquad.x1 = processInputBuffer[i];
    biquad.y2 = biquad.y1;
    biquad.y1 = processOutputBuffer[i];
  } */

  /**
   * Upsample to 48kHz * UPSAMPLE_FACTOR,
   * using polyphase interpolation (79 + 1 zero-pad)-tap FIR filter.
   */
  arm_fir_interpolate_f32(&interp, processOutputBuffer, oversampledOutputBuffer, BUFFER_SIZE);

  /**
   * Scaling buffer back to DAC's 12 bit format, centered at 2048.
   * Values are cast to int32_t for more headroom, clamped, then
   * cast to uint16_t to avoid wraparound distortion
   */
  for (int i = 0; i < UPSAMPLE_BUFFER_SIZE; ++i)
  {
    int32_t out_clamped = (int32_t)(oversampledOutputBuffer[i] * 2048.0f + 2048.0f);
    if (out_clamped > 4095)
      out_clamped = 4095;
    if (out_clamped < 0)
      out_clamped = 0;
    output_buffer[i] = (uint16_t)out_clamped;
  }
  HAL_GPIO_WritePin(GPIOA, GPIO_PIN_0, GPIO_PIN_RESET);
}

/* USER CODE END 4 */

/**
 * @brief  This function is executed in case of error occurrence.
 * @retval None
 */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}
#ifdef USE_FULL_ASSERT
/**
 * @brief  Reports the name of the source file and the source line number
 *         where the assert_param error has occurred.
 * @param  file: pointer to the source file name
 * @param  line: assert_param error line source number
 * @retval None
 */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
