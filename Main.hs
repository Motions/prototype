{-# LANGUAGE RecordWildCards #-}
import Prototype
import Options.Applicative
import Types
import System.IO
import System.Random

data Settings = Settings
    { settingsChainLength :: !Int 
    , settingsLaminsFile :: !FilePath
    , settingsBindersFile :: !FilePath
    , settingsOutputFile :: !FilePath
    , settingsRadius :: !Double
    , settingsNumBinders :: !Int
    , settingsNumSteps :: !Int
    }

parser :: Parser Settings
parser = Settings
    <$> option auto
        (long "chain-length"
        <> short 'l'
        <> help "Chain length"
        <> value 512)
    <*> argument str
        (metavar "LAMIN-BSITES"
        <> help "File containing the lamin binding sites")
    <*> argument str
        (metavar "BINDER-BSITES"
        <> help "File containing the regular binding sites")
    <*> option str
        (long "output"
         <> short 'o'
         <> metavar "OUTPUT"
         <> help "The output file")
    <*> option auto
        (long "radius"
        <> short 'r'
        <> metavar "RADIUS"
        <> help "Radius of the bounding sphere"
        <> value 10)
    <*> option auto
        (long "num-binders"
        <> short 'b'
        <> metavar "BINDERS"
        <> help "Number of binders"
        <> value 256)
    <*> option auto
        (long "steps"
        <> short 's'
        <> metavar "STEPS"
        <> help "Number of simulaton steps"
        <> value 100000)

makeInput :: Settings -> IO Input
makeInput Settings{..} = do
    let inputChainLength = settingsChainLength
        inputRadius = settingsRadius
        inputNumBinders = settingsNumBinders
        inputNumSteps = settingsNumSteps
    inputLamins <- map read . lines <$> readFile settingsLaminsFile
    inputBinders <- map read . lines <$> readFile settingsBindersFile
    inputRandGen <- newStdGen
    return Input{..}

main :: IO ()
main = do
    settings <- execParser $ info (helper <*> parser)
        (fullDesc
        <> progDesc "Perform a MCMC simulation of chromatine movements")
    outputFile <- openFile (settingsOutputFile settings) WriteMode
    input <- makeInput settings
    writePDB outputFile (Prototype.run input)

