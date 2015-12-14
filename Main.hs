{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TemplateHaskell #-}
import Bio.Motions.Prototype as Prototype
import Control.Monad.State.Strict
import Options.Applicative
import System.IO
import System.Random
import Control.Lens

data Settings = Settings
    { settingsChainLength :: !Int 
    , settingsLaminsFile :: !FilePath
    , settingsBindersFile :: !FilePath
    , settingsOutputFile :: !FilePath
    , settingsRadius :: !Int
    , settingsNumBinders :: !Int
    , settingsNumSteps :: !Int
    , settingsWriteIntermediateStates :: !Bool
    }

data Counter = Counter
    { _counterNumSteps :: !Int
    , _counterNumFrames :: !Int
    }
makeLenses ''Counter

parser :: Parser Settings
parser = Settings
    <$> option auto
        (long "chain-length"
        <> short 'l'
        <> help "Chain length"
        <> value 512)
    <*> strArgument
        (metavar "LAMIN-BSITES"
        <> help "File containing the lamin binding sites")
    <*> strArgument
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
    <*> switch
        (long "intermediate-states"
        <> short 'i'
        <> help "Whether to write the intermediate states to the output file")

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

runAndWrite :: (MonadState SimulationState m, MonadIO m) => Maybe Handle -> StateT Counter m ()
runAndWrite handle = do
    counterNumSteps += 1

    oldEnergy <- lift $ gets energy
    lift simulateStep
    newEnergy <- lift $ gets energy

    when (oldEnergy /= newEnergy) $ pushPDB handle

pushPDB :: (MonadState SimulationState m, MonadIO m) => Maybe Handle -> StateT Counter m ()
pushPDB Nothing = return ()
pushPDB (Just handle) = do
    st@SimulationState{..} <- lift get

    headerSequenceNumber <- use counterNumFrames
    headerStep <- use counterNumSteps
    let headerTitle = "chromosome;bonds=" ++ show energy

    liftIO $ do
        putStrLn $ "gyration radius: " ++ show gyrationRadius
        putStrLn $ "energy:          " ++ show energy
        writePDB handle Header{..} st
        hPutStrLn handle "END"

    counterNumFrames += 1

main :: IO ()
main = do
    settings <- execParser $ info (helper <*> parser)
        (fullDesc
        <> progDesc "Perform a MCMC simulation of chromatine movements")
    input <- makeInput settings
    withFile (settingsOutputFile settings) WriteMode $ \outputFile ->
        case initialize input of
            Left e -> print e
            Right st ->
                void $ flip execStateT st $ flip execStateT (Counter 0 0) $ do
                    replicateM_ (settingsNumSteps settings) $ runAndWrite $
                        if settingsWriteIntermediateStates settings then
                            Just outputFile
                        else
                            Nothing
                    pushPDB $ Just outputFile
